#!/usr/bin/env python3
"""Build a native 64-bit candidate in a fresh directory, without installation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=root / "build")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cc", default="gcc")
    parser.add_argument("--domain-offset", action="store_true",
                        help="link the experimental Clipper2 deposition domain backend")
    parser.add_argument("--sanitize", action="store_true",
                        help="build with AddressSanitizer and UndefinedBehaviorSanitizer")
    parser.add_argument("--ordering", choices=["legacy", "amd"], default="legacy",
                        help="AMD requires SuiteSparse AMD development headers and library")
    parser.add_argument("--fe-solver", choices=["legacy", "umfpack"], default="legacy",
                        help="scaled pivoted finite-element LU; requires SuiteSparse UMFPACK")
    parser.add_argument("--geometry-core", type=Path, default=None,
                        help="CMake build directory of suprem4/core: link the geometry core "
                             "(active at run time with SUPREM4_GEOMETRY_CORE=1)")
    args = parser.parse_args()
    args.build_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="native64_", dir=args.build_root.resolve()))
    print(f"Build directory: {work}", flush=True)
    flags = ("-m64 -g -O1 -fcommon -std=gnu89 -Werror=pointer-to-int-cast "
             "-Werror=int-to-pointer-cast -Werror=int-conversion "
             "-Werror=incompatible-pointer-types")
    if args.ordering == "amd":
        flags += " -DSUPREM_USE_AMD"
    if args.fe_solver == "umfpack":
        flags += " -DSUPREM_FE_UMFPACK"
    if args.domain_offset:
        flags += " -DSUPREM_DOMAIN_OFFSET"
    if args.sanitize:
        flags += " -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie -no-pie"
    # A fresh snapshot prevents linking objects built for another architecture
    # or with different flags. Original source-tree objects remain untouched.
    def ignore(directory, names):
        return [name for name in names
                if name.endswith((".o", ".a", ".str"))
                or name.startswith("core")
                or (Path(directory, name).is_file() and name in ("suprem4", "keyread"))]

    shutil.copytree(root / "src", work / "src", ignore=ignore)
    # Generated Makefiles can retain dependencies on 32-bit system headers.
    shutil.copy2(root / "src/Makefile.proto", work / "src/Makefile")
    shutil.copytree(root / "data", work / "data")
    command = ["make", "-j" + str(max(1, args.jobs)), "suprem4",
               "CC=" + args.cc, "COPT=" + flags]
    libraries = ["-lamd"] if args.ordering == "amd" else []
    if args.fe_solver == "umfpack":
        libraries += ["-lumfpack"]
    offset_commands = []
    if args.domain_offset:
        vendor = work / 'src/third_party/clipper2'
        sources = [work / 'src/refine/domain_offset.cpp', work / 'src/refine/oxide_repair.cpp', *sorted((vendor / 'src').glob('*.cpp'))]
        cxxflags = ['-std=c++17', '-g', '-O1', '-I', str(vendor / 'include')]
        if args.sanitize:
            cxxflags += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-fno-pie']
        for source in sources:
            offset_commands.append(['c++', *cxxflags, '-c', str(source), '-o', str(source.with_suffix('.o'))])
        archive = work / 'src/libdomain_offset.a'
        offset_commands.append(['ar', 'rcs', str(archive), *map(lambda p:str(p.with_suffix('.o')), sources)])
        libraries += [str(archive), '-lstdc++']
    if args.geometry_core:
        core = args.geometry_core.resolve()
        # The core (with Mmg) is pre-linked into one object that exports only
        # its C API.  It carries its own Clipper2 / CDT predicates, compiled
        # with other flags than SUPREM's domain-offset copy; linked as plain
        # archives, the linker resolved SUPREM's deposition to the core's copy
        # and a build with the core switched off was no longer bit-identical
        # to one without it (WP7).
        api = ['suprem_core_oxidation_check', 'suprem_core_set_param']
        combined, exported = work / 'src/suprem_core_all.o', work / 'src/suprem_core_api.o'
        offset_commands.append(['ld', '-r', '--force-group-allocation', '-o', str(combined), '--whole-archive',
                                str(core / 'libsuprem_core.a'), str(core / 'third_party/mmg/lib/libmmg2d.a'),
                                '--no-whole-archive'])
        offset_commands.append(['objcopy', *[f'--keep-global-symbol={s}' for s in api],
                                str(combined), str(exported)])
        libraries = [str(exported), *libraries]
        if '-lstdc++' not in libraries:
            libraries.append('-lstdc++')
    command.append('LIBS=' + ' '.join([*libraries, '-lm']))
    manifest = {"command": command, "source_root": str(root),
                "offset_commands": offset_commands,
                "source_sha256": {str(p.relative_to(work / "src")):
                                  hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in (work / "src").rglob("*")
                                  if p.is_file() and (p.suffix in (".c", ".h", ".hpp", ".cpp", ".f", ".proto")
                                                     or p.name in ("Makefile", "LICENSE"))},
                "data_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (work / "data").iterdir() if p.is_file()}}
    (work / "build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log = work / "build.log"
    with log.open("w") as output:
        for offset_command in offset_commands:
            result = subprocess.run(offset_command, cwd=work / 'src', stdout=output, stderr=subprocess.STDOUT)
            if result.returncode:
                print(f'Domain backend build failed: {log}')
                return 1
        result = subprocess.run(command, cwd=work / "src", stdout=output,
                                stderr=subprocess.STDOUT)
    if result.returncode:
        print(f"Build failed; diagnostics: {log}")
        print("\n".join(log.read_text().splitlines()[-35:]))
        return 1
    binary = work / "src/suprem4"
    ident = binary.read_bytes()[:5]
    if ident != b"\x7fELF\x02":
        print(f"Expected a 64-bit ELF executable: {binary}")
        return 1
    print(f"Native 64-bit candidate: {binary}\nBuild log: {log}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
