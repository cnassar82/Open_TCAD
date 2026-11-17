/* Real binary key reader: first card, nested siblings, and return to root. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "shell.h"
#include "check.h"
struct command_table command[NUMCMD];
int get_proc(char *name) { return 0; }
extern int read_ukfile(char *);
static void record(FILE *fd, int level, char *name, int value)
{
    char fixed[12] = {0};
    int type = INT, zero = 0;
    strncpy(fixed, name, 11);
    assert(fwrite(&level, sizeof(int), 1, fd) == 1);
    fwrite(fixed, sizeof(fixed), 1, fd);
    fwrite(&type, sizeof(int), 1, fd);
    fwrite(&value, sizeof(int), 1, fd);
    fputc(0, fd); fputc(0, fd);
    fwrite(&zero, sizeof(int), 1, fd);
}
static void release(struct par_str **list)
{
    int i;
    for (i=0; list[i]; ++i) {
        if (list[i]->param) { release(list[i]->param); free(list[i]->param); }
        free(list[i]); list[i] = NULL;
    }
}
int main(int argc, char **argv)
{
    FILE *fd;
    assert(argc == 2);
    fd = fopen(argv[1], "wb"); assert(fd);
    record(fd, 0, "first", 10);
    record(fd, 1, "child", 20);
    record(fd, 2, "grandchild", 30);
    record(fd, 1, "sibling", 40);
    record(fd, 0, "second", 50);
    fclose(fd);
    unsetenv("KEYFILE");
    assert(read_ukfile(argv[1]) == 1);
    assert(!strcmp(cards[0]->name, "first"));
    assert(cards[0]->def.ival == 10);
    assert(cards[0]->param[0]->def.ival == 20);
    assert(cards[0]->param[0]->param[0]->def.ival == 30);
    assert(cards[0]->param[0]->param[1] == NULL);
    assert(cards[0]->param[1]->def.ival == 40);
    assert(cards[0]->param[2] == NULL);
    assert(cards[1]->def.ival == 50 && cards[2] == NULL);
    release(cards);
    puts("PASS key reader first/nested/sibling/root records");
    return 0;
}
