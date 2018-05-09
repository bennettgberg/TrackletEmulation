#include "tracklet_em_2.h"
/*
 * Version 2018.2: March 29, 2018
 *
 */

etaphibin * L1_cluster(etaphibin * phislice){

		etaphibin * clusters = (etaphibin *)malloc(netabins/2 * sizeof(etaphibin));
	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		float my_pt, left_pt, right_pt, right2pt, left2pt;
		int nclust = 0;

		for(int etabin = 0; etabin < netabins; ++etabin){
			//assign values for my pT and neighbors' pT
			my_pt = phislice[etabin].pTtot;
			if(etabin > 0 && !phislice[etabin-1].used) {
				left_pt = phislice[etabin-1].pTtot;
				if(etabin > 1 && !phislice[etabin-2].used) {
					left2pt = phislice[etabin-2].pTtot;
				} else left2pt = 0;
			} else left_pt = 0;
			if(etabin < netabins - 1 && !phislice[etabin+1].used) {
				right_pt = phislice[etabin+1].pTtot;
				if(etabin < netabins - 2 && !phislice[etabin+2].used) {
					right2pt = phislice[etabin+2].pTtot;
				} else right2pt = 0;
			} else right_pt = 0;
		
		//if I'm not a cluster, move on.
			if(my_pt < left_pt || my_pt <= right_pt) {
			   //if unused pT in the left neighbor, spit it out as a cluster.
			        if(left_pt > 0) {
					clusters[nclust] = phislice[etabin-1];
					phislice[etabin-1].used = true;
					++nclust;
				}
				continue;
			}

		//I guess I'm a cluster-- should I use my right neighbor?
		// Note: left neighbor will definitely be used because if it 
		//       didn't belong to me it would have been used already
			clusters[nclust] = phislice[etabin];
			phislice[etabin].used = true;
			if(left_pt > 0) {
				clusters[nclust].pTtot += left_pt;
				clusters[nclust].numtracks += phislice[etabin-1].numtracks;
			}
			if(my_pt >= right2pt && right_pt > 0) {
				clusters[nclust].pTtot += right_pt;
				clusters[nclust].numtracks += phislice[etabin+1].numtracks;
				phislice[etabin+1].used = true;
			}

			++nclust;
		} //for each etabin                       
	                         
	//Now merge clusters, if necessary
		for(int m = 0; m < nclust -1; ++m){
			if(fabs(clusters[m+1].eta - clusters[m].eta) < 1.5*etastep){
				if(clusters[m+1].pTtot > clusters[m].pTtot){
					clusters[m].eta = clusters[m+1].eta;
				}
				clusters[m].pTtot += clusters[m+1].pTtot;
				for(int m1 = m+1; m1 < nclust-1; ++m1){
					clusters[m1] = clusters[m1+1];
				}
				nclust--;
				m = -1;
			}//end if clusters neighbor in eta
		}//end for (m) loop
//	for(int i = 0; i < nclust; ++i) cout << clusters[i].phi << "\t" << clusters[i].pTtot << "\t" << clusters[i].numtracks << endl;
	//zero out remaining unused clusters.
	for(int i = nclust; i < netabins/2; ++i){
		clusters[i].pTtot = 0;
	}
	return clusters;
}
