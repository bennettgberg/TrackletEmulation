#include "tracklet_em_2.h"

//  for 1trk, file called "trk_output_v2.txt" should contain all event data in format:
//"Event" [2 spaces] # [3 spaces] jet pT [3 spaces] jet eta [3 spaces] jet phi  (for both jets)
// pT [3 spaces] eta [3 spaces] phi [3 spaces] z     (for each track)
//
// For example: 
// 
// Event  3   173.579330444   0.999113619328   2.01476287842
// Event  3   173.579330444   -0.999113619328   -1.12682986259
// 19.2072582245   0.990153074265   2.02103352547   -3.41288733482
// ...
// Then, next event.
//
//  for tp, file called "tp_output_v2.txt" in same format
//Arguments: argv[1]: 1 for 1trk, 2 for tp. argv[2]: number of zbins (default 16).
int main(int argc, char ** argv){
   //Numbers to specify which events to start and finish on.
	int eventstart = 0;
	int eventend = 99999;      //if higher than number of events, will end with last event.
   
   //If second number specified it is number of zbins
   //Default 16.
	string nz = "16";
	int nzbins = 16;
	if(argc == 3){
		nz = argv[2];
		nzbins = atoi(argv[2]);
	}
	else if(argc > 3) {
		cout << "Error: too many arguments" << endl;
		exit(0);
	}
	else if(argc < 2){
		cout << "Error: please specify '1' for 1trk or '2' for tp." << endl;
		exit(0);
	}
    //Name of input file (1st argument: 1 is 1trk, 2 is tp).	
	string filename;
        if(atoi(argv[1]) == 1){
 	    filename = "trk_output_v2.txt";
        }
        else if(atoi(argv[1]) == 2){
            filename = "tp_output_v2.txt";
        } 
        else {
            cout << "Error: please specify 1 for 1trk or 2 for tp."<<endl;
            exit(0);
        }

	track_data * tracks = (track_data *)malloc(numtracks * sizeof(track_data));
	string data_in;
  //Open input file and all output files.
	ifstream in_tracks;
	in_tracks.open(filename.c_str());
	ofstream out_clusts;
	string outname;
        if(atoi(argv[1]) == 1){
             outname = "em_out_trk_" + nz + "z.txt";
        } 
        else {
             outname = "em_out_tp_" + nz + "z.txt";
        }
	out_clusts.open(outname.c_str());
	ofstream phidat;
	string plotname = "phidat_" + nz + "z.dat";
	phidat.open(plotname.c_str());
	ofstream etadat;
	plotname = "etadat_" + nz + "z.dat";
	etadat.open(plotname.c_str());
	ofstream ptdat;
	plotname = "ptplot_" + nz + "z.dat";
	ptdat.open(plotname.c_str());
	ofstream traxdat;
	plotname = "ntracks_" + nz + "z.dat";
	traxdat.open(plotname.c_str());
	ofstream disdat;
	plotname = "distance_" + nz + "z.dat";
	disdat.open(plotname.c_str());
	ofstream zplot;
	plotname = "zplot_" + nz + "z.dat";
	zplot.open(plotname.c_str());
        ofstream phitie;
        phitie.open("phi_pt.dat");
//	ofstream trackphi;
  //      trackphi.open("trackphi.dat");
	int nevents = 0;
	getline(in_tracks, data_in, ' ');
        string data;
	float distance;
     	while(nevents <= eventend){
 //pointer to structure to hold jet MC data for each event.
        	mc_data * mcdat = (mc_data *)malloc(10*sizeof(mc_data));
		int ntracks = 0;
		int ntp = 0;
		bool use_event = true;
		string evnum = "0";
		while(data_in == "Event"){
			for(int cc=0; cc<2; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			if(data_in != evnum){
				++nevents;
				ntp = 0;
			}
			evnum = data_in;
			for(int cc=0; cc<2; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			getline(in_tracks, data, ' ');
                        mcdat[ntp].ogpt = atof(data.c_str());
			for(int cc=0; cc<3; ++cc){
				getline(in_tracks, data, ' ');
			}
                        mcdat[ntp].ogeta = atof(data.c_str());
			getline(in_tracks, data);
                        data = data.substr(0, data.length()-1);
                        mcdat[ntp].ogphi = atof(data.c_str());
			mcdat[ntp].ntracks = 0;                //start at 0, find total later.
			ntp++;
		      //Now get next particle (or first track data)
			getline(in_tracks, data_in, ' ');
		}
		if(nevents < eventstart){
			ntracks = 0;
			while(data_in != "Event"){
				getline(in_tracks, data_in, ' ');
			}
			continue;
		}
	   	out_clusts << "\n****EVENT " << nevents << " *** (" << ntp << " tracking particles)" << endl;
		for(int i = 0; i < ntp; ++i){
	   		out_clusts << i << " pT: " << mcdat[i].ogpt << " eta: " << mcdat[i].ogeta << " phi: " << mcdat[i].ogphi << endl;
		}
		//getline(in_tracks, data_in, ' ');
		if(in_tracks.eof()){
			//cout << "End of file reached!" << endl;
			break;
		}
		while(data_in != "Event") {
		//Convert strings to floating point numbers!!
			tracks[ntracks].pT = atof(data_in.c_str());
			//read spaces
			for(int cc=0; cc<3; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			tracks[ntracks].eta = atof(data_in.c_str());
			for(int cc=0; cc<3; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			tracks[ntracks].phi = atof(data_in.c_str());
			//read spaces
			for(int cc=0; cc<2; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			getline(in_tracks, data_in);
			tracks[ntracks].z = atof(data_in.c_str());
			tracks[ntracks].bincount = 0;
			if(tracks[ntracks].pT >= 2.0 && fabs(tracks[ntracks].eta) < 2.4 && fabs(tracks[ntracks].z) < 15.0){
				ntracks++;
			}    
			else {
				use_event = false;
			}
		//Read next line of data
			getline(in_tracks, data_in, ' ');
			if(in_tracks.eof()){
				break;
			}

                } //data_in isn't "Event"
		if(!use_event){
			out_clusts << "EVENT DISCARDED (DATA OUT OF RANGE)" << endl;
			continue;
		}
  //              for(int t=0; t < ntracks; ++t){
//			trackphi << tracks[t].phi << "\t" << tracks[t].pT << endl;
//
//		}

		//do matching for input tracks to tracking particles
		int bbb = 0;
		for(int b = 0; b < ntracks; b++){
			float mindist = 15;
	     		for(int k = 0; k < ntp; ++k){
		   		distance = sqrt(pow(mcdat[k].ogphi - tracks[b].phi, 2) + pow(mcdat[k].ogeta - tracks[b].eta, 2));
				if(distance < mindist){
					bbb = k;
					mindist = distance;
				}
			}//for each cluster in the zbin
			//cout << "Min distance for track " << b << ": " << mindist << endl;
			mcdat[bbb].ntracks++;
		}//for each track


          //find clusters for tracks data.
             out_clusts << ntracks << " tracks total" << endl;
	//     cout << "Event " << nevents << ": " << ntp << " tracking particles, " << ntracks << " tracks " << endl;
             maxzbin * mzb = L2_cluster(tracks, mcdat, nzbins, ntracks);
             if(mzb == NULL){
                continue;
             }
         //output to plot data files.
           for(int zz=0; zz < nzbins-1; ++zz){
             out_clusts << "ZBIN: " << all_zbins[zz].znum << endl; 
	     if(&all_zbins[zz] == mzb) out_clusts << "VERTEX" << endl;
             for(int k = 0; k < all_zbins[zz].nclust; ++k){
          //     cout << all_zbins[zz].clusters[k].phi <<"\t" << all_zbins[zz].clusters[k].pTtot << "\t" << all_zbins[zz].clusters[k].numtracks << endl;//print data on tying phibins together.
		float mindist = 15.0;
                int closeboi = -1;
		bool matched = false;
         	for(int b = 0; b < ntp; b++){
               //Match clusters with correct input jet (and ignore the garbage clusters.)
               // find minimum distance from cluster to track.
               // Only accept cluster if it is within .3 (in eta-phi space) of one of the jets.
                
		   distance = sqrt(pow(all_zbins[zz].mcd[b].ogphi - all_zbins[zz].clusters[k].phi, 2) + pow(all_zbins[zz].mcd[b].ogeta - all_zbins[zz].clusters[k].eta, 2));
		   if(distance < mindist) {
			mindist = distance;
			closeboi = b;
		   }
		} //for each input particle
		   if (closeboi >= 0 && mindist < 0.3){
	     		phidat << all_zbins[zz].mcd[closeboi].ogphi << "\t" << all_zbins[zz].clusters[k].phi << endl;
			etadat << all_zbins[zz].mcd[closeboi].ogeta << "\t" << all_zbins[zz].clusters[k].eta << endl;
                 //in pt data file: print phi, MC pT, cluster pT, then number of tracks of the cluster.
			ptdat << all_zbins[zz].clusters[k].phi << "\t" << all_zbins[zz].mcd[closeboi].ogpt << "\t" << all_zbins[zz].clusters[k].pTtot << "\t" << all_zbins[zz].clusters[k].numtracks << endl;
			traxdat << ntracks << "\t" <<  all_zbins[zz].clusters[k].numtracks << endl;
			disdat << all_zbins[zz].mcd[closeboi].ogpt << "\t" << distance << endl;
			out_clusts << "(match dist=" << mindist <<") CLUSTER " << k << "(" << all_zbins[zz].clusters[k].numtracks << " tracks)\t MC data \t Matched cluster data" << endl;
			out_clusts << "   phi \t\t\t\t\t\t" << all_zbins[zz].mcd[closeboi].ogphi << "\t\t" << all_zbins[zz].clusters[k].phi << endl;
			out_clusts << "   eta \t\t\t\t\t\t" << all_zbins[zz].mcd[closeboi].ogeta << "\t\t" << all_zbins[zz].clusters[k].eta << endl;
			out_clusts << "   pT \t\t\t\t\t\t" << all_zbins[zz].mcd[closeboi].ogpt << "\t\t" << all_zbins[zz].clusters[k].pTtot << endl;
			matched = true;
               cout << all_zbins[zz].mcd[closeboi].ogphi <<"\t" << all_zbins[zz].mcd[closeboi].ogpt << "\t" <<all_zbins[zz].clusters[k].pTtot << "\t" << all_zbins[zz].clusters[k].numtracks << endl;//print data on tying phibins together.
			continue; //break;
	     	   } 
	        if(matched)
			continue;
		out_clusts << "  (Unmatched) CLUSTER " << k << "(" << all_zbins[zz].clusters[k].numtracks << " tracks)" << endl;
		out_clusts << "   phi \t\t\t\t\t\t\t\t" << all_zbins[zz].clusters[k].phi << endl;
		out_clusts << "   eta \t\t\t\t\t\t\t\t" << all_zbins[zz].clusters[k].eta << endl;
		out_clusts << "   pT \t\t\t\t\t\t\t\t"  << all_zbins[zz].clusters[k].pTtot << endl;
	     } //for each cluster
            } //for each zbin
	    
		
	     for(int b = 0; b < ntp; b++){
	     	float mindist = 15;
		etaphibin match = mzb->clusters[0];
		for(int zb = 0; zb < nzbins-1; zb++){
			if(&all_zbins[zb] == NULL){
				continue;
			}
			for(int k = 0; k < all_zbins[zb].nclust; ++k){
		   		distance = sqrt(pow(mcdat[b].ogphi - all_zbins[zb].clusters[k].phi, 2) + pow(mcdat[b].ogeta - all_zbins[zb].clusters[k].eta, 2));
				if(distance < mindist){
					mindist = distance;
					match = all_zbins[zb].clusters[k];
				}
				
			}//for each cluster in the zbin
		}//for each zbin
		//print data to output files ***********************************	
		zplot << mcdat[b].ntracks << "\t" << match.numtracks << "\t" << mcdat[b].ogpt << "\t" << match.pTtot << endl;
	     }
	
	for(int i = 0; i < nzbins-1; i++){	
	  if(&all_zbins[i] == NULL || &all_zbins[i].mcd == NULL) continue;
	    free(all_zbins[i].clusters);
	}	  
   	free(all_zbins);
	free(mcdat);

	} //while nevents <= eventend
        out_clusts << "****" << nevents << " events****" << endl;
	free(tracks);
	zplot.close();
        in_tracks.close();
	out_clusts.close();
	phidat.close();
	etadat.close();
	ptdat.close();
	traxdat.close();
	disdat.close();
        phitie.close();
//        trackphi.close();
//	cout << "6 new data files created." <<endl;
	return 0;
} //end main
