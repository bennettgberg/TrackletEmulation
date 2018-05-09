# TrackletEmulation
Updated 11/10/2017: now separates cluster-finding into two layers to more accurately emulate the firmware.
L1_cluster does phi-clustering in each etabin, and L2_cluster stitches everything together.

11/21/2017: updated to use an algorithm closer to the firmware algorithm. Also some bug fixes.

11/22/2017: fixed bug with counting number of tracks

1/16/2018: updated L2_cluster algorithm to include jet-merging at end, some minor bug fixes

3/29/2018: Updated to do clustering first in eta (L1_cluster), then in phi (L2_cluster).
           This includes connecting first phibin to the last one. Possibly still some bugs (updates to come).
     
4/2/2018: Simplified algorithm to tie two sides of phi together (to reduce the logic in the firmware)

5/9/2018: Added 3 versions of L2 clustering: no tying of first+last phibins, simple tying, and complex tying.
          Also some bug fixes.
