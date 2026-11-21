#ifndef THIN_FILM_H
#define THIN_FILM_H
/* Enforce at least `layers` element layers across every dielectric film
 * (dbase/thin_film.c). Returns the number of conforming edge bisections. */
extern int thin_film_conform(int layers);
#endif
