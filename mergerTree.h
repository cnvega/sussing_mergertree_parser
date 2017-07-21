#ifndef MERGERTREE_H
#define MERGERTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#include <inttypes.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#define LINEMAX 1024

struct Halo {
   int snapshot;
   int64_t id;
   int nProgenitors;
   int64_t firstProgenitor;
   int64_t sibling;
};

struct Halo **HaloSnap;
int64_t     *HaloStoreIdx;

FILE *Forests, *MergerTree;

int64_t TotNumHalos;
int64_t *NumHalosSnap;
int     MaxSnapshot;

/* It obtains the total number of halos in each snapshot */
void forest_stats(char *, int);

/* It opens the forest list and merger tree files, and it allocates the needed
 * memory to store the complete MT structure. */
void open_catalogs(char *, int);

/* Process the next forest by reading it and assigning the historic links of
 * the halos. It returns the ID of the processed forests or zero if it has 
 * finished. */
int64_t forest_process();

/* It reorders all the loaded halos according to their ids (in each
 * snapshot) */
void reorder_by_ids();

/* It saves the complete parsed merger tree in a hdf5 file. */
herr_t save_mergertree(char *);
#endif
