#ifndef MERGERTREE_C
#define MERGERTREE_C

#include "mergerTree.h"

/* Local prototypes */
static int order_by_ids(const void *, const void *);
static struct Halo* find_halo_by_id(struct Halo*, int64_t, int64_t);

/* Public functions */

void open_catalogs(char *path, int maxSnapshot)
{
   char fname[LINEMAX];
   char line[LINEMAX], ctmp;
   fpos_t fpos;
   FILE *stats;
   int snap;
   int64_t i64tmp;
  
    
   sprintf(fname, "nhalos.list");
   if(!(stats = fopen(fname, "r"))) {
      fprintf(stderr, "ERROR: you need to create this file first: `%s`\n", fname);
      exit(EXIT_FAILURE);
   }
   fscanf(stats, "%"SCNd64, &TotNumHalos);
   NumHalosSnap = malloc((maxSnapshot+1)*sizeof(int64_t));
   for (snap=0; snap<=maxSnapshot; snap++)
      fscanf(stats, "%"SCNd64, &NumHalosSnap[snap]);
   fclose(stats);


   sprintf(fname, "%s/sussing_forests.list", path);
   if(!(Forests = fopen(fname, "r"))) {
      fprintf(stderr, "ERROR: can't open file `%s`\n", fname);
      exit(EXIT_FAILURE);
   }
   // skipping the lines with comments:
   fgetpos(Forests, &fpos);
   ctmp = fgetc(Forests);
   while (ctmp == '#')
   {
      fgets(line, LINEMAX, Forests);
      fgetpos(Forests, &fpos);
      ctmp = fgetc(Forests);
   }
   fsetpos(Forests, &fpos);

   
   sprintf(fname, "%s/sussing_tree.list", path);
   
   if(!(MergerTree = fopen(fname, "r")))
   {
      fprintf(stderr, "ERROR: can't open file `%s`\n", fname);
      exit(EXIT_FAILURE);
   }
   while (fgetc(MergerTree) != '\n');  // Format version
   while (fgetc(MergerTree) != '\n');  // Description
   // Check the total number of halos
   fscanf(MergerTree, "%"SCNd64, &i64tmp);
   
   if (i64tmp != TotNumHalos) {
      fprintf(stderr, "ERROR: The file nhalos.list is not consistent with: `%s`\n", fname);
      exit(EXIT_FAILURE);
   }

   // Now we can allocate the memory:
   HaloSnap = malloc((maxSnapshot+1)*sizeof(struct Halo*));
   HaloStoreIdx = malloc((maxSnapshot+1)*sizeof(int64_t));
   
   for(snap=0; snap<=maxSnapshot; snap++) {
      HaloStoreIdx[snap] = 0;
      HaloSnap[snap] = calloc(NumHalosSnap[snap], sizeof(struct Halo));
   }
   MaxSnapshot = maxSnapshot;

}

int64_t forest_process()
{
   int64_t fid, nhalos;
   int64_t *nharr;
   int snap, start;
  
   int64_t hid, progId, k;
   struct Halo *firstp, *lastp;
   int     nprog, p;

   struct Halo **halosForest, h={0};

   nharr = malloc((MaxSnapshot+1)*sizeof(int64_t));
   halosForest = malloc((MaxSnapshot+1)*sizeof(struct Halo*));

   // load and process the stuff
   if (fscanf(Forests, "%"SCNd64, &fid) == EOF)
      return 0;
   
   fscanf(Forests, "%"SCNd64, &nhalos);
   
   for (snap=0; snap<=MaxSnapshot; snap++) {
      fscanf(Forests, "%"SCNd64, &nharr[snap]);
      halosForest[snap] = malloc(nharr[snap]*sizeof(struct Halo));
   }
   fscanf(Forests, "\n");
   
   for (snap=0; snap<=MaxSnapshot; snap++) {
      for (k=0; k<nharr[snap]; k++) {
         fscanf(MergerTree, "%"SCNd64" %"SCNd64"\n", &h.id, &h.nProgenitors);
         h.snapshot = snap;
         h.sibling = h.id;
         if (h.nProgenitors > 0) {            
            // reading the first progenitor first:
            fscanf(MergerTree, "%"SCNd64"\n", &h.firstProgenitor);
            firstp = find_halo_by_id(&halosForest[snap-1][0], nharr[snap-1], h.firstProgenitor);
            lastp = firstp;
            
            // now the other progenitors and create linked list:
            for (p=1; p<h.nProgenitors; p++) {
               fscanf(MergerTree, "%"SCNd64"\n", &progId);
               lastp->sibling = progId;
               lastp = find_halo_by_id(&halosForest[snap-1][0], nharr[snap-1], progId);
            }
            lastp->sibling = h.firstProgenitor;
         }
         else
            h.firstProgenitor = 0;

         halosForest[snap][k] = h;
      }

      // Finally, we need to order this new array:
      if (nharr[snap] > 1)
         qsort(&halosForest[snap][0], nharr[snap], sizeof(struct Halo), order_by_ids);
   }
   
   // save them to the corresponding pointers.
   for (snap=0; snap<=MaxSnapshot; snap++) { 
      for (k=0; k<nharr[snap]; k++) {
         HaloSnap[snap][HaloStoreIdx[snap]++] = halosForest[snap][k];
      }
      free(halosForest[snap]);
   }
   free(halosForest);

   
   return fid;
}

void reorder_by_ids()
{
   int snap;
   // Order all the arrays at the end:
   for (snap=0; snap<=MaxSnapshot; snap++)
      qsort(&HaloSnap[snap][0], NumHalosSnap[snap], sizeof(struct Halo), order_by_ids);
}

herr_t save_mergertree(char *path)
{
   char fname[LINEMAX], label[LINEMAX], grlabel[LINEMAX];
   int snap;

   int64_t *i64buf, i;
   int     *ibuf;

   hid_t file_id, group_id;
   herr_t status;
   hsize_t dims[2] = {0, 1};

   sprintf(fname, "%s/sussing_tree_ids.h5", path);
   
   file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

   // create the group tree:
   for (snap=0; snap<=MaxSnapshot; snap++) {
      sprintf(grlabel, "snapshot_%03d", snap);
      group_id = H5Gcreate(file_id, grlabel, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Gclose(group_id);

      dims[0] = NumHalosSnap[snap];
      
      i64buf = malloc(NumHalosSnap[snap]*sizeof(int64_t));
      ibuf = malloc(NumHalosSnap[snap]*sizeof(int));

      sprintf(label, "%s/HaloID", grlabel);
      for (i=0; i<NumHalosSnap[snap]; i++)
         i64buf[i] = HaloSnap[snap][i].id;
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT64, i64buf);
      
      sprintf(label, "%s/NumProgenitors", grlabel);
      for (i=0; i<NumHalosSnap[snap]; i++)
         ibuf[i] = HaloSnap[snap][i].nProgenitors;
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT, ibuf);
      
      sprintf(label, "%s/FirstProgenitor", grlabel);
      for (i=0; i<NumHalosSnap[snap]; i++)
         i64buf[i] = HaloSnap[snap][i].firstProgenitor;
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT64, i64buf);
      
      sprintf(label, "%s/Sibling", grlabel);
      for (i=0; i<NumHalosSnap[snap]; i++)
         i64buf[i] = HaloSnap[snap][i].sibling;
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_INT64, i64buf);

      free(i64buf);
      free(ibuf);
   }
   status = H5Fclose(file_id);

   printf("File %s saved.\n", fname);
}

void forest_stats(char *path, int maxSnapshot)
{
   char fname[LINEMAX];
   char ctmp, line[LINEMAX];
   fpos_t fpos;
   FILE *fin, *fout;
   int snap;

   int64_t id, nhalos, *nharr;
   int64_t nhalosTotal=0, *nharrTotal;
   
   nharr = malloc((maxSnapshot+1)*sizeof(int64_t));
   nharrTotal = malloc((maxSnapshot+1)*sizeof(int64_t));

   for (snap=0; snap<=maxSnapshot; snap++)
      nharrTotal[snap] = 0;

   sprintf(fname, "%s/sussing_forests.list", path);
   
   if(!(fin = fopen(fname, "r")))
   {
      fprintf(stderr, "ERROR: can't open file `%s`\n", fname);
      exit(EXIT_FAILURE);
   }
   /* skip lines with comments: */
   fgetpos(fin, &fpos);
   ctmp = fgetc(fin);
   while (ctmp == '#')
   {
      fgets(line, LINEMAX, fin);
      fgetpos(fin, &fpos);
      ctmp = fgetc(fin);
   }
   fsetpos(fin, &fpos);

   while(fscanf(fin, "%"SCNd64, &id) != EOF)
   {
      fscanf(fin, "%"SCNd64, &nhalos);
      
      nhalosTotal += nhalos;
      //printf("%"PRId64" %"PRId64" ", id, nhalos);
      for (snap=0; snap<=maxSnapshot; snap++)
      {
         fscanf(fin, "%"SCNd64, &nharr[snap]);
         nharrTotal[snap] += nharr[snap];
         //printf("%"PRId64" ", nharr[snap]);
      }
      //printf("\n");
      fscanf(fin, "\n");
   }

   fclose(fin);
   
   fout = fopen("nhalos.list", "w");
   fprintf(fout, "%"PRId64"\n", nhalosTotal);
   for (snap=0; snap<=maxSnapshot; snap++)
      fprintf(fout, "%"PRId64"\n", nharrTotal[snap]);
   fclose(fout);

}


/* Local functions: */

// Sorting criterion: 
int order_by_ids(const void *a, const void *b)
{
   const struct Halo *aa = a;
   const struct Halo *bb = b;
   if (aa->id < bb->id) return -1;
   if (aa->id > bb->id) return 1;
   return 0;
}

// Search in array:
struct Halo* find_halo_by_id(struct Halo *hlist, int64_t nhalos, int64_t hid)
{
   int64_t i;
   int64_t min=0, max=nhalos-1;
   if (hid < 0) return 0;
   while (min < max)
   {
      i = min + (max-min)/2;
      if (hlist[i].id < hid) min = i+1;
      else max = i;
   }
   if (hlist[max].id == hid) return &hlist[max];

   // return -1;
   fprintf(stderr, "ERROR: Id %"PRId64" not found in array!\n\n", hid);
   exit(EXIT_FAILURE);
}

#endif
