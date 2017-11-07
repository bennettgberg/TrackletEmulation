#include "tracklet_em.h"

//input array of track_data, output zbin of maximum ht.
struct maxzbin * find_clusters(struct track_data * tracks, struct mc_data * mcd, int nzbins){
    //returns NULL if there are no tracks for this event.
        int ntracks = mcd->ntracks;
        if(ntracks == 0){
	      return NULL;
	}
        const float zstep = 2.0 * maxz / nzbins;
        
	float zmin = -1.0*maxz;
	float zmax = zmin + 2*zstep;
	//Create grid of etaphibins! 
	etaphibin epbins[nphibins][netabins];
	float phi = -1.0 * M_PI;
	float eta;
	for(int i = 0; i < nphibins; ++i){
	    eta = -1.0 * maxeta;
            for(int j = 0; j < netabins; ++j){
		epbins[i][j].phimin = phi;
		epbins[i][j].phimax = phi + phistep;
		epbins[i][j].etamin = eta;
		eta = eta + etastep;
		epbins[i][j].etamax = eta;
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
	        for(int i = 0; i < nphibins; ++i){
			for(int j = 0; j < netabins; ++j){
				epbins[i][j].pTtot = 0;
				epbins[i][j].used = false;
				epbins[i][j].numtracks = 0;
			}//for each etabin
		} //for each phibin

	      //Fill in etaphibins grid with pT from each track.
		for(int k = 0; k < ntracks; ++k) {
			for(int i = 0; i < nphibins; ++i){
				for(int j = 0; j < netabins; ++j){
					if((zmin <= tracks[k].z && zmax >= tracks[k].z) &&
					  ((epbins[i][j].etamin <= tracks[k].eta && epbins[i][j].etamax >= tracks[k].eta) 
					    && epbins[i][j].phimin <= tracks[k].phi && epbins[i][j].phimax >= tracks[k].phi && tracks[k].bincount != 2)){
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
	

		//Create clusters array to hold output cluster data; can't have more clusters than tracks.
		struct etaphibin * clusters = (struct etaphibin *)malloc(ntracks * sizeof(struct etaphibin));

	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		int hipT;
		int nmax;
		int pmax;
		int nclust = 0;

		while(true){
		 	hipT = 0;
			//Find eta-phibin with highest pT.
			for(int etabin = 0; etabin < netabins; ++etabin){
				for(int phibin = 0; phibin < nphibins; ++phibin){
					if(!epbins[phibin][etabin].used && epbins[phibin][etabin].pTtot >= hipT){
						hipT = epbins[phibin][etabin].pTtot;
						nmax = etabin;
						pmax = phibin;
					}
				}//for each phibin
			}//for each etabin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			clusters[nclust] = epbins[pmax][nmax];
			epbins[pmax][nmax].used = true;

		      //Add pT of the 8 neighbors.
			if(nmax + 1 < netabins && !epbins[pmax][nmax + 1].used){
				clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax][nmax + 1].pTtot;
				epbins[pmax][nmax + 1].used = true;
				clusters[nclust].numtracks += epbins[pmax][nmax+1].numtracks;
			}

			if(nmax - 1 >= 0 && !epbins[pmax][nmax - 1].used){
				clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax][nmax - 1].pTtot;
				epbins[pmax][nmax - 1].used = true;
				clusters[nclust].numtracks += epbins[pmax][nmax-1].numtracks;
			}

			if(pmax-1 >= 0){
				if(nmax - 1 >= 0 && !epbins[pmax-1][nmax - 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax-1][nmax - 1].pTtot;
					epbins[pmax-1][nmax - 1].used = true;
					clusters[nclust].numtracks += epbins[pmax-1][nmax-1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(nmax - 1 >= 0 && !epbins[nphibins-1][nmax - 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[nphibins-1][nmax - 1].pTtot;
					epbins[nphibins-1][nmax - 1].used = true;
					clusters[nclust].numtracks += epbins[nphibins-1][nmax-1].numtracks;
				}
	
			if(pmax-1 >= 0){
				if(nmax - 1 >= 0 && !epbins[pmax-1][nmax].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax-1][nmax].pTtot;
					epbins[pmax-1][nmax].used = true;
					clusters[nclust].numtracks += epbins[pmax-1][nmax].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(nmax - 1 >= 0 && !epbins[nphibins-1][nmax].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[nphibins-1][nmax].pTtot;
					epbins[nphibins-1][nmax].used = true;
					clusters[nclust].numtracks += epbins[nphibins-1][nmax].numtracks;
				}

			if(pmax-1 >= 0){
				if(nmax + 1 < netabins && !epbins[pmax-1][nmax + 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax-1][nmax + 1].pTtot;
					epbins[pmax-1][nmax + 1].used = true;
					clusters[nclust].numtracks += epbins[pmax-1][nmax+1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(nmax + 1 < netabins && !epbins[nphibins-1][nmax + 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[nphibins-1][nmax + 1].pTtot;
					epbins[nphibins-1][nmax + 1].used = true;
					clusters[nclust].numtracks += epbins[nphibins-1][nmax+1].numtracks;
				}

			if(pmax+1 < nphibins){
				if(nmax - 1 >= 0 && !epbins[pmax+1][nmax - 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax+1][nmax - 1].pTtot;
					epbins[pmax+1][nmax - 1].used = true;
					clusters[nclust].numtracks += epbins[pmax+1][nmax-1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(nmax - 1 >= 0 && !epbins[0][nmax - 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[0][nmax - 1].pTtot;
					epbins[0][nmax - 1].used = true;
					clusters[nclust].numtracks += epbins[0][nmax-1].numtracks;
				}

			if(pmax+1 < nphibins){
				if(!epbins[pmax+1][nmax].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax+1][nmax].pTtot;
					epbins[pmax+1][nmax].used = true;
					clusters[nclust].numtracks += epbins[pmax+1][nmax].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(!epbins[0][nmax].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[0][nmax].pTtot;
					epbins[0][nmax].used = true;
					clusters[nclust].numtracks += epbins[0][nmax].numtracks;
				}

			if(pmax+1 < nphibins){
				if(nmax + 1 < netabins && !epbins[pmax+1][nmax + 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[pmax+1][nmax + 1].pTtot;
					epbins[pmax+1][nmax + 1].used = true;
					clusters[nclust].numtracks += epbins[pmax+1][nmax+1].numtracks;
				}
			}
			else  //pmax is on the edge; its neighbor is on the other side.
				if(nmax + 1 < netabins && !epbins[0][nmax + 1].used){
					clusters[nclust].pTtot = clusters[nclust].pTtot + epbins[0][nmax + 1].pTtot;
					epbins[0][nmax + 1].used = true;
					clusters[nclust].numtracks += epbins[0][nmax+1].numtracks;
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
			ht += clusters[k].pTtot;
                }

	   //if ht is larger than previous max, this is the new vertex zbin.
		if(ht > mzb->ht){
			mzb->znum = zbin;
                      //reinitialize clusters array.
                        free(mzb->clusters);
			mzb->clusters = (struct etaphibin *)malloc(nclust*sizeof(struct etaphibin));
			mzb->nclust = nclust;
			for(int k = 0; k < nclust; ++k){
				mzb->clusters[k].phi = (clusters[k].phimin + clusters[k].phimax)/2.0;
				mzb->clusters[k].eta = (clusters[k].etamin + clusters[k].etamax)/2.0;
				mzb->clusters[k].pTtot = clusters[k].pTtot;
				mzb->clusters[k].numtracks = clusters[k].numtracks;
			}
			mzb->ht = ht;
		}
	       //Prepare for next zbin!
		zmin = zmin + zstep;
		zmax = zmax + zstep;
	     } //for each zbin
         
       return mzb;
}
