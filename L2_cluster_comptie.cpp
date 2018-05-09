#include "tracklet_em_2.h"

maxzbin * all_zbins;
//input array of track_data, output zbin of maximum ht.
maxzbin * L2_cluster(track_data * tracks, mc_data * mcd, int nzbins, int ntracks){
    //returns NULL if there are no tracks for this event.
        if(ntracks == 0){
	      return NULL;
	}
	all_zbins = (maxzbin *)malloc(nzbins*sizeof(maxzbin));
        const float zstep = 2.0 * maxz / nzbins;
        
	float zmin = -1.0*maxz;
	float zmax = zmin + 2*zstep;
	//Create grid of phibins! 
	etaphibin ** epbins = (etaphibin **)malloc(nphibins * sizeof(etaphibin *));
	for(int i = 0; i < nphibins; ++i){
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
		epbins[i][j].phi = (phimin + phimax) / 2;
		epbins[i][j].eta = (etamin + etamax) / 2;
	     }//for each phibin
	     phi = phi + phistep;
	 } //for each etabin (finished creating epbins)

	maxzbin * mzb = &all_zbins[0];
	//mzb->mcd = mcd;
        //mzb->clusters = (etaphibin*)malloc(sizeof(etaphibin));
	 //Last zbin won't be used (goes beyond maximum z)
	for(int zbin = 0; zbin < nzbins-1; ++zbin){
	
	      //First initialize pT, numtracks, used to 0 (or false)
	        for(int i = 0; i < nphibins; ++i){
			for(int j = 0; j < netabins; ++j){
				epbins[i][j].pTtot = 0;
				epbins[i][j].used = false;
				epbins[i][j].numtracks = 0;
			}//for each phibin
		} //for each phibin

	      //Fill in etaphibins grid with pT from each track.
		for(int k = 0; k < ntracks; ++k) {
			for(int i = 0; i < nphibins; ++i){
				for(int j = 0; j < netabins; ++j){
					if((zmin <= tracks[k].z && zmax >= tracks[k].z) &&
					  ((epbins[i][j].eta - etastep / 2 <= tracks[k].eta && epbins[i][j].eta + etastep / 2 >= tracks[k].eta) 
					    && epbins[i][j].phi - phistep / 2 <= tracks[k].phi && epbins[i][j].phi + phistep / 2 >= tracks[k].phi && tracks[k].bincount != 2)){
						++tracks[k].bincount;
						epbins[i][j].pTtot += tracks[k].pT;
						++epbins[i][j].numtracks;
					} //if right bin
				} //for each phibin: j loop
			}//for each phibin: i loop
		} //for each track: k loop

    //Uncomment to print out pT of each eta and phi bin.
//		for(int i = 0; i < nphibins; ++i)
//			for(int j = 0; j < nphibins; ++j)
//				if(epbins[i][j].pTtot != 0) {
//				        cout << "zmin " << zmin << " zmax " << zmax << endl;
//					cout << "zbin " << zbin << " epbins[" << i << "][" << j << "] pTtot: " << epbins[i][j].pTtot << endl;
//				}
	

	  //First do clustering in Layer 1: maximum possible nclust for each eta slice would be a cluster in every other phibin.
	   float totalpt = 0;
		etaphibin ** L1clusters = (etaphibin**)malloc(nphibins*sizeof(etaphibin*));
                for(int phislice = 0; phislice < nphibins; ++phislice){
			L1clusters[phislice] = L1_cluster(epbins[phislice]);
			for(int ind = 0; L1clusters[phislice][ind].pTtot != 0; ++ind){
				L1clusters[phislice][ind].used = false;
				totalpt += L1clusters[phislice][ind].pTtot;
			}
		}

	//Create clusters array to hold output cluster data for Layer2; can't have more clusters than tracks.
		etaphibin * L2cluster = (etaphibin *)malloc(ntracks * sizeof(etaphibin));

	//Find eta-phibin with maxpT, make center of cluster, add neighbors if not already used.
		float hipT = 0;
		int nclust = 0;
		int phibin = 0;
		int imax;
	     //index of clusters array for each phislice.
		int index1;
		float E1 =0;
		float E0 =0;
		float E2 =0;
		int trx1, trx2;
		int used1, used2, used3, used4;

		bool last_used[netabins];
		bool first_used[netabins];
		bool pre_seen[netabins], pre_seen2[netabins];
		for(int eb = 0; eb < netabins; eb++) {
			last_used[eb] = false;
			first_used[eb] = false;
			pre_seen[eb] = false;
			pre_seen2[eb] = false;
		}
		/* Tie together first and last phibin! */
		//cluster centered at phibin 0?
		while(true){
			hipT = 0;
			for(index1 = 0; L1clusters[nphibins-1][index1].pTtot > 0; ++index1){
				if(!pre_seen[index1] && !L1clusters[nphibins-1][index1].used && L1clusters[nphibins-1][index1].pTtot >= hipT){
					hipT = L1clusters[nphibins-1][index1].pTtot;
					imax = index1;
				}
			}//for each index within the phibin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			E0 = hipT;
			E1 = 0;
			E2 = 0;
			trx1 = 0;
			trx2 = 0;
			used1 = -1;
			used2 = -1;
		    //now find pT of phibin 0, store in E1.
			for(int im = 0; L1clusters[0][im].pTtot != 0; ++im){
				if(L1clusters[0][im].used == true) continue;
				if(fabs(L1clusters[0][im].eta - L1clusters[0][imax].eta) < 1.5*etastep){
					E1 += L1clusters[0][im].pTtot;
					trx1 += L1clusters[0][im].numtracks;
					if(used1 >= 0) used2 = im; else used1 = im;
				}
			}				
			pre_seen[imax] = true;
			if(E1 == 0) continue;   //just added!
		// if last <= phibin 0, cluster @ phibin 0, use last.
			if(E0 <= E1) {
			//	L2cluster[nclust] = L1clusters[nphibins-1][imax];
				L2cluster[nclust].phi = L1clusters[0][used1].phi;
				L2cluster[nclust].pTtot = E0+E1;
				L2cluster[nclust].numtracks += trx1;
				L1clusters[0][used1].used = true;
				if(used2 >= 0) L1clusters[0][used2].used = true;
				L1clusters[nphibins-1][imax].used = true;
			}
			used3 = -1;
			used4 = -1;
                //now get pT of phibin 1, store in E2.
			for(int im = 0; L1clusters[1][im].pTtot != 0; ++im){
				if(L1clusters[1][im].used == true) continue;
				if(fabs(L1clusters[1][im].eta - L1clusters[1][imax].eta) < 1.5*etastep){
					E2 += L1clusters[1][im].pTtot;
					trx2 += L1clusters[1][im].numtracks;
					if(used3 >= 0) used4 = im; else used3 = im;
				}
			}
		//if the cluster is centered at phibin 0, using last, just see if we should add E2 or nah.
			if(E0 <= E1) {
				if(E2 < E1) {
					L2cluster[nclust].pTtot += E2;
					L2cluster[nclust].numtracks += trx2;
					L1clusters[1][used3].used = true;
					if(used4 >= 0) L1clusters[1][used4].used = true;
				}
				++nclust;
				continue;
			}
	//if phibin 1 >= last phibin, there is a cluster centered either at phibin 0 or phibin 1 (but not using last)
			if(E2 >= E0) {
				int trx0 = 0;
				int used5 = -1; int used6 = -1;
                          //if phibin 0 > phibin 1, then phibin 0 is center of the cluster.
				if(E1 > E2) {
					L2cluster[nclust] = L1clusters[0][used1];
				}
				else {
					L2cluster[nclust] = L1clusters[1][used3];
				    //get pT of second phibin, see if it should be added.
				    //  now E0 will be for phibin 2.
					E0 = 0;
					for(int im = 0; L1clusters[2][im].pTtot != 0; ++im){
						if(L1clusters[2][im].used == true) continue;
						if(fabs(L1clusters[2][im].eta - L1clusters[2][imax].eta) < 1.5*etastep){
							E2 += L1clusters[2][im].pTtot;
							trx0 += L1clusters[2][im].numtracks;
							if(used5 >= 0) used6 = im; else used5 = im;
						}
					}
				//if we're not using phibin 2, set E0 and trx0 to 0.
					if(E0 >= E1){
						 E0 = 0;
						 trx0 = 0;
					}
				}
				L2cluster[nclust].eta = L1clusters[nphibins-1][imax].eta;
				L2cluster[nclust].pTtot = E1 + E2 + E0; //E0 is only non-zero if phibin 1 < phibin 0.
				L2cluster[nclust].numtracks = trx1 + trx2 + trx0;
				L1clusters[0][used1].used = true;
				if(used2 >= 0) L1clusters[0][used2].used = true;
				L1clusters[1][used3].used = true;
				if(used4 >= 0) L1clusters[1][used4].used = true;
				if(used5 >= 0) L1clusters[2][used5].used = true;
				if(used6 >= 0) L1clusters[2][used6].used = true;
				++nclust;
			}
		}
                	

		/* Now finish clustering the rest of the phibins. */

			//Find eta-phibin with highest pT.
		for(phibin = 0; phibin < nphibins; ++phibin){
		    while(true){
			hipT = 0;
			for(index1 = 0; L1clusters[phibin][index1].pTtot > 0; ++index1){
				if(!L1clusters[phibin][index1].used && L1clusters[phibin][index1].pTtot >= hipT){
					hipT = L1clusters[phibin][index1].pTtot;
					imax = index1;
				}
			}//for each index within the phibin
		      //If highest pT is 0, all bins are used.
			if(hipT == 0){
				break;
			}
			E0 = hipT;
			E1 = 0;
			E2 = 0;
			trx1 = 0;
			trx2 = 0;
			L2cluster[nclust] = L1clusters[phibin][imax];
			L1clusters[phibin][imax].used = true;
		//Add pT of neighbors.
		//Higher neighbors(s)
			if(phibin != nphibins-1){
				used1 = -1;
				used2 = -1;
				for (index1 = 0; L1clusters[phibin+1][index1].pTtot != 0; ++index1){
					if(L1clusters[phibin+1][index1].used){
						continue;
					}
					if(fabs(L1clusters[phibin+1][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
						E1 += L1clusters[phibin+1][index1].pTtot;
						trx1 += L1clusters[phibin+1][index1].numtracks;
						if(used1 < 0)
							used1 = index1;
						else
							used2 = index1;
					}//if cluster is within one phibin
		
				} //end for each cluster in above phibin
				if(E1 < E0){
					L2cluster[nclust].pTtot += E1;   
					L2cluster[nclust].numtracks += trx1;
					if(used1 >= 0)
						L1clusters[phibin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[phibin+1][used2].used = true;
					++nclust;
					continue;
				}
				
				if(phibin != nphibins-2){
					used3 = -1;
					used4 = -1;
				
					for (index1 = 0; L1clusters[phibin+2][index1].pTtot != 0; ++index1){
						if(L1clusters[phibin+2][index1].used){
							continue;
						}
						if(fabs(L1clusters[phibin+2][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
							E2 += L1clusters[phibin+2][index1].pTtot;
							trx2 += L1clusters[phibin+2][index1].numtracks;
							if(used3 < 0)
								used3 = index1;
							else
								used4 = index1;
						}
		
					}
					if(E2 < E1){
						L2cluster[nclust].pTtot += E1 + E2;
						L2cluster[nclust].numtracks += trx1 + trx2;
						L2cluster[nclust].eta = L1clusters[phibin+1][used1].eta;	
						if(used1 >= 0)
							L1clusters[phibin+1][used1].used = true;
						if(used2 >= 0)
							L1clusters[phibin+1][used2].used = true;
						if(used3 >= 0)
							L1clusters[phibin+2][used3].used = true;
						if(used4 >= 0)
							L1clusters[phibin+2][used4].used = true;
					}
					++nclust;
				}
				else { //if it is second-to-last phibin (25)
				//find pT of phibin 0
					used3 = -1;
					used4 = -1;
					for (index1 = 0; L1clusters[0][index1].pTtot != 0; ++index1){
						if(L1clusters[0][index1].used){
							continue;
						}
						if(fabs(L1clusters[0][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
							E2 += L1clusters[0][index1].pTtot;
							trx2 += L1clusters[0][index1].numtracks;
							if(used3 < 0)
								used3 = index1;
							else
								used4 = index1;
						}
					}
				//add pT of all 3 phibins together
					L2cluster[nclust].pTtot += E1 + E2;
					L2cluster[nclust].numtracks += trx1 + trx2;
					L2cluster[nclust].eta = L1clusters[phibin+1][used1].eta;	
					if(used1 >= 0)
						L1clusters[phibin+1][used1].used = true;
					if(used2 >= 0)
						L1clusters[phibin+1][used2].used = true;
					if(used3 >= 0)
						L1clusters[0][used3].used = true;
					if(used4 >= 0)
						L1clusters[0][used4].used = true;
					++nclust;
				}
			}//end Not last phibin(26)
			else { //if it is last phibin (26)
				//L1clusters[phibin][imax].used = true;
				//find pT of first phibin (if unused) 
				used1 = -1;
				used2 = -1;
				for (index1 = 0; L1clusters[0][index1].pTtot != 0; ++index1){
					if(L1clusters[0][index1].used){
						continue;
					}
					if(fabs(L1clusters[0][index1].eta - L1clusters[phibin][imax].eta) <= 1.5*etastep){
						E1 += L1clusters[0][index1].pTtot;
						trx1 += L1clusters[0][index1].numtracks;
						if(used1 < 0)
							used1 = index1;
						else
							used2 = index1;
					}//if cluster is within one phibin
				} //for each cluster in last phibin
			//whatever pT is left in phibin 0 definitely belongs to last.
				L2cluster[nclust].pTtot += E1;   
				L2cluster[nclust].numtracks += trx1;
				if(used1 >= 0)
					L1clusters[0][used1].used = true;
				if(used2 >= 0)
					L1clusters[0][used2].used = true;
				++nclust;
			}
		    }//while hipT not 0
		}//for each phibin
	
	//Now merge clusters, if necessary
		for(int m = 0; m < nclust -1; ++m){
                     for(int n = m+1; n < nclust; ++n)
                        if(L2cluster[n].eta == L2cluster[m].eta && (fabs(L2cluster[n].phi - L2cluster[m].phi) < 1.5*phistep /*|| fabs(L2cluster[n].phi - L2cluster[m].phi) > 6.0*/)){
                                if(L2cluster[n].pTtot > L2cluster[m].pTtot){
                                        L2cluster[m].eta = L2cluster[n].eta;
                                }
                                L2cluster[m].pTtot += L2cluster[n].pTtot;
                                L2cluster[m].numtracks += L2cluster[n].numtracks;
                                for(int m1 = n; m1 < nclust-1; ++m1){
                                        L2cluster[m1] = L2cluster[m1+1];
                                }
                                nclust--;
                                m = -1;
                                break; //?????
                        }//end if clusters neighbor in eta
                }//end for (m) loop     
//		for(int db=0;db<nclust;++db)cout<<L2cluster[db].phi<<"\t"<<L2cluster[db].pTtot<<"\t"<<L2cluster[db].numtracks<<endl;	
          //sum up all pTs in this zbin to find ht.
		float ht = 0;
		for(int k = 0; k < nclust; ++k){
			ht += L2cluster[k].pTtot;
                }

	   //if ht is larger than previous max, this is the new vertex zbin.
	   	all_zbins[zbin].mcd = mcd;
		all_zbins[zbin].znum = zbin;
		all_zbins[zbin].clusters = (etaphibin *)malloc(nclust*sizeof(etaphibin));
		all_zbins[zbin].nclust = nclust;
		for(int k = 0; k < nclust; ++k){
			all_zbins[zbin].clusters[k].phi = L2cluster[k].phi;                               
			all_zbins[zbin].clusters[k].eta = L2cluster[k].eta;                             
			all_zbins[zbin].clusters[k].pTtot = L2cluster[k].pTtot;
			all_zbins[zbin].clusters[k].numtracks = L2cluster[k].numtracks;
		}
		all_zbins[zbin].ht = ht;
		if(ht > mzb->ht){
			mzb = &all_zbins[zbin];
		}
	       //Prepare for next zbin!
		zmin = zmin + zstep;
		zmax = zmax + zstep;
	     } //for each zbin
       return mzb;
}
