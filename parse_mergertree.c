#include "mergerTree.h"

int main()
{
   char path[LINEMAX];
   int maxSnapshot;

   int64_t fid=0;
   
   maxSnapshot = 125;
   sprintf(path, "/foo/bar/sussing/simulation/path");
   
   open_catalogs(path, maxSnapshot);

   while(forest_process()) {
      fid++;
   };
   
   printf("Computed %"PRId64" forests\n", fid);

   reorder_by_ids();
   save_mergertree(".");

   return EXIT_SUCCESS;
}

