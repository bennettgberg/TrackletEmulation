#include "tracklet_em.h"

struct etaphibin * L1_cluster(struct etaphibin * etaslice){

		struct etaphibin * clusters = (struct etaphibin *)malloc(nphibins/2 * sizeof(struct etaphibin));
	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		int hipT;
		int pmax;
		int nclust = 0;

		while(true){
		 	hipT = 0;
			//Find eta-phibin with highest pT.
			for(int phibin = 0; phibin < nphibins; ++phibin){
				if(!etaslice[phibin].used && etaslice[phibin].pTtot >= hipT){
					hipT = etaslice[phibin].pTtot;
					pmax = phibin;
				}
			}//for each phibin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			clusters[nclust] = etaslice[pmax];
			etaslice[pmax].used = true;

		      //Add pT of the 2 neighbors.
		      //Lower neighbor
			if(pmax-1 >= 0){
				if(!etaslice[pmax-1].used){
					clusters[nclust].pTtot += etaslice[pmax-1].pTtot;
					etaslice[pmax-1].used = true;
					clusters[nclust].numtracks += etaslice[pmax-1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(!etaslice[nphibins-1].used){
					clusters[nclust].pTtot += etaslice[nphibins-1].pTtot;
					etaslice[nphibins-1].used = true;
					clusters[nclust].numtracks += etaslice[nphibins-1].numtracks;
				}
                            	
                      //Higher neighbor
			if(pmax+1 < nphibins){
				if(!etaslice[pmax+1].used){
					clusters[nclust].pTtot += etaslice[pmax+1].pTtot;
					etaslice[pmax+1].used = true;
					clusters[nclust].numtracks += etaslice[pmax+1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(!etaslice[0].used){
					clusters[nclust].pTtot += etaslice[0].pTtot;
					etaslice[0].used = true;
					clusters[nclust].numtracks += etaslice[0].numtracks;
				}

			++nclust;
		} //while hipT isn't 0 (still unused bins)
	//zero out remaining unused clusters.
	for(int i = nclust; i < nphibins/2; ++i){
		clusters[i].pTtot = 0;
	}
	return clusters;
}
