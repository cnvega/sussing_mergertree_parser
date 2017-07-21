#include "mergerTree.h"

int main()
{
   char path[LINEMAX];
   int maxSnapshot;
   
   maxSnapshot = 125;
   sprintf(path, "/foo/bar/sussing/simulation/path");
   
   forest_stats(path, maxSnapshot);

   return EXIT_SUCCESS;
}
