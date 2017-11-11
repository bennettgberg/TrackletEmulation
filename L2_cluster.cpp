#include "tracklet_em.h"

//input array of track_data, output zbin of maximum ht.
struct maxzbin * L2_cluster(struct track_data * tracks, struct mc_data * mcd, int nzbins){
    //returns NULL if there are no tracks for this event.
        int ntracks = mcd->ntracks;
        if(ntracks == 0){
	      return NULL;
	}
        const float zstep = 2.0 * maxz / nzbins;
        
	float zmin = -1.0*maxz;
	float zmax = zmin + 2*zstep;
	//Create grid of phibins! 
	struct etaphibin ** epbins = (struct etaphibin **)malloc(netabins * sizeof(struct etaphibin *));
	for(int i = 0; i < netabins; ++i){
		epbins[i] = (struct etaphibin *)malloc(nphibins * sizeof(struct etaphibin));
	}
	float phi = -1.0 * M_PI;
	float eta;
	float etamin, etamax, phimin, phimax;
	for(int i = 0; i < nphibins; ++i){
	    eta = -1.0 * maxeta;
            for(int j = 0; j < netabins; ++j){
		phimin = phi;
		phimax = phi + phistep;
		etamin = eta;
		eta = eta + etastep;
		etamax = eta;
		epbins[j][i].phi = (phimin + phimax) / 2;
		epbins[j][i].eta = (etamin + etamax) / 2;
	     }//for each etabin
	     phi = phi + phistep;
	 } //for each phibin (finished creating epbins)

	maxzbin * mzb = new maxzbin();
	mzb->phimc = mcd->ogphi;         
	mzb->etamc = mcd->ogeta;          
	mzb->pTmc = mcd->ogpt;
        mzb->clusters = (struct etaphibin*)malloc(sizeof(etaphibin));
	 //Last zbin won't be used (goes beyond maximum z)
	for(int zbin = 0; zbin < nzbins-1; ++zbin){
	
	      //First initialize pT, numtracks, used to 0 (or false)
	        for(int i = 0; i < netabins; ++i){
			for(int j = 0; j < nphibins; ++j){
				epbins[i][j].pTtot = 0;
				epbins[i][j].used = false;
				epbins[i][j].numtracks = 0;
			}//for each etabin
		} //for each phibin

	      //Fill in etaphibins grid with pT from each track.
		for(int k = 0; k < ntracks; ++k) {
			for(int i = 0; i < netabins; ++i){
				for(int j = 0; j < nphibins; ++j){
					if((zmin <= tracks[k].z && zmax >= tracks[k].z) &&
					  ((epbins[i][j].eta - etastep / 2 <= tracks[k].eta && epbins[i][j].eta + etastep / 2 >= tracks[k].eta) 
					    && epbins[i][j].phi - phistep / 2 <= tracks[k].phi && epbins[i][j].phi + phistep / 2 >= tracks[k].phi && tracks[k].bincount != 2)){
						++tracks[k].bincount;
						epbins[i][j].pTtot += tracks[k].pT;
						++epbins[i][j].numtracks;
					} //if right bin
				} //for each etabin: j loop
			}//for each phibin: i loop
		} //for each track: k loop

    //Uncomment to print out pT of each eta and phi bin.
	//	for(int i = 0; i < nphibins; ++i)
	//		for(int j = 0; j < netabins; ++j)
	//			cout << "epbins[" << i << "][" << j << "] pTtot: " << epbins[i][j].pTtot << endl;
	

	  //First do clustering in Layer 1: maximum possible nclust for each eta slice would be a cluster in every other phibin.
		struct etaphibin ** L1clusters = (struct etaphibin**)malloc(netabins*sizeof(struct etaphibin*));
                for(int etaslice = 0; etaslice < netabins; ++etaslice){
			L1clusters[etaslice] = L1_cluster(epbins[etaslice]);
			for(int ind = 0; L1clusters[etaslice][ind].pTtot != 0; ++ind){
				L1clusters[etaslice][ind].used = false;
			}
		}

	//Create clusters array to hold output cluster data for Layer2; can't have more clusters than tracks.
		struct etaphibin * L2clusters = (struct etaphibin *)malloc(ntracks * sizeof(struct etaphibin));

	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		float hipT = 0;
		int nmax;
		int nclust = 0;
		int etabin = 0;
		int imax;
	     //index of clusters array for each etaslice.
		int index1;
		while(true){
		 	hipT = 0;
			//Find eta-phibin with highest pT.
			for(etabin = 0; etabin < netabins; ++etabin){
				for(index1 = 0; L1clusters[etabin][index1].pTtot > 0; ++index1){
					if(!L1clusters[etabin][index1].used && L1clusters[etabin][index1].pTtot >= hipT){
						hipT = L1clusters[etabin][index1].pTtot;
						nmax = etabin;
						imax = index1;
					}
				}//for each index within the etabin
			}//for each etabin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			L2clusters[nclust] = L1clusters[nmax][imax];
			L1clusters[nmax][imax].used = true;
		//Add pT of neighbors.
		//Lower neighbor(s)
			if(nmax != 0){
				for (index1 = 0; L1clusters[nmax-1][index1].pTtot != 0; ++index1){
					if(L1clusters[nmax-1][index1].used){
						continue;
					}
					if(fabs(L1clusters[nmax-1][index1].phi - L1clusters[nmax][imax].phi) <= phistep){
						L2clusters[nclust].pTtot += L1clusters[nmax][imax].pTtot;
						L1clusters[nmax-1][index1].used = true;
					}
				}

			}

		//Higher neighbors(s)
			if(nmax != netabins-1){
				for (index1 = 0; L1clusters[nmax+1][index1].pTtot != 0; ++index1){
					if(L1clusters[nmax+1][index1].used){
						continue;
					}
					if(fabs(L1clusters[nmax+1][index1].phi - L1clusters[nmax][imax].phi) <= phistep){
						L2clusters[nclust].pTtot += L1clusters[nmax][imax].pTtot;
						L1clusters[nmax+1][index1].used = true;
					}
				}

			}

			++nclust;
		} //while hipT isn't 0 (still unused bins)

	      //if no clusters in this zbin, nothing to print
		if(nclust == 0){
			zmin = zmin + zstep;
			zmax = zmax + zstep;
			continue;
		}
	
          //sum up all pTs in this zbin to find ht.
		float ht = 0;
		for(int k = 0; k < nclust; ++k){
			ht += L2clusters[k].pTtot;
                }

	   //if ht is larger than previous max, this is the new vertex zbin.
		if(ht > mzb->ht){
			mzb->znum = zbin;
                      //reinitialize clusters array.
                        free(mzb->clusters);
			mzb->clusters = (struct etaphibin *)malloc(nclust*sizeof(struct etaphibin));
			mzb->nclust = nclust;
			for(int k = 0; k < nclust; ++k){
				mzb->clusters[k].phi = L2clusters[k].phi;                               
				mzb->clusters[k].eta = L2clusters[k].eta;                             
				mzb->clusters[k].pTtot = L2clusters[k].pTtot;
				mzb->clusters[k].numtracks = L2clusters[k].numtracks;
			}
			mzb->ht = ht;
		}
	       //Prepare for next zbin!
		zmin = zmin + zstep;
		zmax = zmax + zstep;
	     } //for each zbin
         
       return mzb;
}
