#include "tracklet_em_2.h"

//input array of track_data, output zbin of maximum ht.
maxzbin * L2_cluster(track_data * tracks, mc_data * mcd, int nzbins){
    //returns NULL if there are no tracks for this event.
        int ntracks = mcd->ntracks;
        if(ntracks == 0){
	      return NULL;
	}
        const float zstep = 2.0 * maxz / nzbins;
        
	float zmin = -1.0*maxz;
	float zmax = zmin + 2*zstep;
	//Create grid of phibins! 
	etaphibin ** epbins = (etaphibin **)malloc(netabins * sizeof(etaphibin *));
	for(int i = 0; i < netabins; ++i){
		epbins[i] = (etaphibin *)malloc(nphibins * sizeof(etaphibin));
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
/*	mzb->phimc = mcd->ogphi;         
	mzb->etamc = mcd->ogeta;          
	mzb->pTmc = mcd->ogpt;*/
	mzb->mcd = mcd;
        mzb->clusters = (etaphibin*)malloc(sizeof(etaphibin));
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
//		for(int i = 0; i < netabins; ++i)
//			for(int j = 0; j < nphibins; ++j)
//				if(epbins[i][j].pTtot != 0) {
//				        cout << "zmin " << zmin << " zmax " << zmax << endl;
//					cout << "zbin " << zbin << " epbins[" << i << "][" << j << "] pTtot: " << epbins[i][j].pTtot << endl;
//				}
	

	  //First do clustering in Layer 1: maximum possible nclust for each eta slice would be a cluster in every other phibin.
		etaphibin ** L1clusters = (etaphibin**)malloc(netabins*sizeof(etaphibin*));
                for(int etaslice = 0; etaslice < netabins; ++etaslice){
			L1clusters[etaslice] = L1_cluster(epbins[etaslice]);
			for(int ind = 0; L1clusters[etaslice][ind].pTtot != 0; ++ind){
				L1clusters[etaslice][ind].used = false;
			}
		}

	//Create clusters array to hold output cluster data for Layer2; can't have more clusters than tracks.
		etaphibin * L2cluster = (etaphibin *)malloc(ntracks * sizeof(etaphibin));

	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		float hipT = 0;
		int nclust = 0;
		int etabin = 0;
		int imax;
	     //index of clusters array for each etaslice.
		int index1;
		float E1 =0;
		float E0 =0;
		float E2 =0;
		int trx1, trx2;
		int used1, used2, used3, used4;

			//Find eta-phibin with highest pT.
		for(etabin = 0; etabin < netabins; ++etabin){
		    while(true){
			hipT = 0;
			for(index1 = 0; L1clusters[etabin][index1].pTtot > 0; ++index1){
				if(!L1clusters[etabin][index1].used && L1clusters[etabin][index1].pTtot >= hipT){
					hipT = L1clusters[etabin][index1].pTtot;
					imax = index1;
				}
			}//for each index within the etabin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			E0 = hipT;
			E1 = 0;
			E2 = 0;
			trx1 = 0;
			trx2 = 0;
			L2cluster[nclust] = L1clusters[etabin][imax];
			L1clusters[etabin][imax].used = true;
		//Add pT of neighbors.
		//Higher neighbors(s)
			if(etabin != netabins-1){
				used1 = -1;
				used2 = -1;
				for (index1 = 0; L1clusters[etabin+1][index1].pTtot != 0; ++index1){
					if(L1clusters[etabin+1][index1].used){
						continue;
					}
					if(fabs(L1clusters[etabin+1][index1].phi - L1clusters[etabin][imax].phi) <= phistep){
						E1 += L1clusters[etabin+1][index1].pTtot;
						trx1 += L1clusters[etabin+1][index1].numtracks;
						if(used1 < 0)
							used1 = index1;
						else
							used2 = index1;
					}//if cluster is within one phibin
		
				} //for each cluster in above etabin
				if(E1 < E0){
					L2cluster[nclust].pTtot += E1;   
					L2cluster[nclust].numtracks += trx1;
					if(used1 >= 0)
						L1clusters[etabin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[etabin+1][used2].used = true;
					++nclust;
					continue;
				}
				
				if(etabin != netabins-2){
					used3 = -1;
					used4 = -1;
					for (index1 = 0; L1clusters[etabin+2][index1].pTtot != 0; ++index1){
						if(L1clusters[etabin+2][index1].used){
							continue;
					}
						if(fabs(L1clusters[etabin+2][index1].phi - L1clusters[etabin][imax].phi) <= phistep){
							E2 += L1clusters[etabin+2][index1].pTtot;
							trx2 += L1clusters[etabin+2][index1].numtracks;
							if(used3 < 0)
								used3 = index1;
							else
								used4 = index1;
						}
		
					}
					if(E2 < E1){
						L2cluster[nclust].pTtot += E1 + E2;
						L2cluster[nclust].numtracks += trx1 + trx2;
						L2cluster[nclust].eta = L1clusters[etabin+1][used1].eta;	
						if(used1 >= 0)
							L1clusters[etabin+1][used1].used = true;
						if(used2 >= 0)
							L1clusters[etabin+1][used2].used = true;
						if(used3 >= 0)
							L1clusters[etabin+2][used3].used = true;
						if(used4 >= 0)
							L1clusters[etabin+2][used4].used = true;
					}
					++nclust;
					continue;
				}
				else{
					L2cluster[nclust].pTtot += E1;
					L2cluster[nclust].numtracks += trx1;
					L2cluster[nclust].eta = L1clusters[etabin+1][used1].eta;
					if(used1 >= 0)
						L1clusters[etabin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[etabin+1][used2].used = true;
					++nclust;
					continue;
				}
			}//Not last etabin
		    }//while hipT not 0
		}//for each etabin
	      //if no clusters in this zbin, nothing to print
		if(nclust == 0){
			zmin = zmin + zstep;
			zmax = zmax + zstep;
			continue;
		}
	
          //sum up all pTs in this zbin to find ht.
		float ht = 0;
		for(int k = 0; k < nclust; ++k){
			ht += L2cluster[k].pTtot;
                }

	   //if ht is larger than previous max, this is the new vertex zbin.
		if(ht > mzb->ht){
			mzb->znum = zbin;
                      //reinitialize clusters array.
                        free(mzb->clusters);
			mzb->clusters = (etaphibin *)malloc(nclust*sizeof(etaphibin));
			mzb->nclust = nclust;
			for(int k = 0; k < nclust; ++k){
				mzb->clusters[k].phi = L2cluster[k].phi;                               
				mzb->clusters[k].eta = L2cluster[k].eta;                             
				mzb->clusters[k].pTtot = L2cluster[k].pTtot;
				mzb->clusters[k].numtracks = L2cluster[k].numtracks;
			}
			mzb->ht = ht;
		}
	       //Prepare for next zbin!
		zmin = zmin + zstep;
		zmax = zmax + zstep;
	     } //for each zbin
       return mzb;
}
