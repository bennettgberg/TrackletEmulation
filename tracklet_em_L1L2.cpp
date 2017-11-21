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
//Arguments: argv[1]: 1 for 1trk, 2 for tp. argv[2]: number of zbins (default 64).
int main(int argc, char ** argv){
   //Numbers to specify which events to start and finish on.
	int eventstart = 0;
	int eventend = 99999;      //if higher than number of events, will end with last event.
   
   //If second number specified it is number of zbins
   //Default 64.
	string nz = "64";
	int nzbins = 64;
	if(argc == 3){
		nz = argv[2];
		nzbins = atoi(argv[2]);
	}
	else if(argc > 3) {
		cout << "Error: too many arguments" << endl;
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

	struct track_data * tracks = (struct track_data *)malloc(numtracks * sizeof(struct track_data));
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
	ofstream plot1dat;
	string plotname = "plot1dat_" + nz + "z.dat";
	plot1dat.open(plotname.c_str());
	ofstream plot2dat;
	plotname = "plot2dat_" + nz + "z.dat";
	plot2dat.open(plotname.c_str());
	ofstream ptdat;
	plotname = "ptplot_" + nz + "z.dat";
	ptdat.open(plotname.c_str());
	ofstream traxdat;
	plotname = "ntracks_" + nz + "z.dat";
	traxdat.open(plotname.c_str());
	ofstream disdat;
	plotname = "distance_" + nz + "z.dat";
	disdat.open(plotname.c_str());
	int nevents = 0;
	getline(in_tracks, data_in, ' ');
        string data;
	float distance;
     	while(nevents < eventend){
 //pointer to structure to hold jet MC data for each event.
        	struct mc_data * mcdat = (struct mc_data *)malloc(10*sizeof(struct mc_data));
		int ntracks = 0;
		int ntp = 0;
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
			cout << "End of file reached!" << endl;
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
		//Read next line of data
			if(in_tracks.eof()){
				cout << "End of file reached!" << endl;
				break;
			}
			getline(in_tracks, data_in, ' ');

             } //data_in isn't "Event"

          //find clusters for tracks data.
             out_clusts << ntracks << " tracks total" << endl;
             mcdat->ntracks = ntracks;
             struct maxzbin * mzb = L2_cluster(tracks, mcdat, nzbins);
             if(mzb == NULL){
                continue;
             }
         //output to plot data files.
             out_clusts << "ZBIN: " << mzb->znum << endl; 
             for(int k = 0; k < mzb->nclust; ++k){
		bool matched = false;
         	for(int b = 0; b < ntp; b++){
               //Match clusters with correct input jet (and ignore the garbage clusters.)
               // Only accept cluster if it is within .3 (in eta-phi space) of one of the jets.
                
		   distance = sqrt(pow(mzb->mcd[b].ogphi - mzb->clusters[k].phi, 2) + pow(mzb->mcd[b].ogeta - mzb->clusters[k].eta, 2));
		   if (distance < 0.3){
	     		plot1dat << mzb->mcd[b].ogphi << "\t" << mzb->clusters[k].phi << endl;
			plot2dat << mzb->mcd[b].ogeta << "\t" << mzb->clusters[k].eta << endl;
			ptdat << mzb->mcd[b].ogpt << "\t" << mzb->clusters[k].pTtot << endl;
			traxdat << ntracks << "\t" <<  mzb->clusters[k].numtracks << endl;
			disdat << mzb->mcd[b].ogpt << "\t" << distance << endl;
			out_clusts << "   CLUSTER " << k << "(" << mzb->clusters[k].numtracks << " tracks)\t MC data \t Matched cluster data" << endl;
			out_clusts << "   phi \t\t\t" << mzb->mcd[b].ogphi << "\t\t" << mzb->clusters[k].phi << endl;
			out_clusts << "   eta \t\t\t" << mzb->mcd[b].ogeta << "\t\t" << mzb->clusters[k].eta << endl;
			out_clusts << "   pT \t\t\t" << mzb->mcd[b].ogpt << "\t\t" << mzb->clusters[k].pTtot << endl;
			matched = true;
			break;
	     	   } 
                }//for each input jet
	        if(matched)
			continue;
		out_clusts << "  (Unmatched) CLUSTER " << k << "(" << mzb->clusters[k].numtracks << " tracks)" << endl;
		out_clusts << "   phi \t\t\t" << "\t\t" << mzb->clusters[k].phi << endl;
		out_clusts << "   eta \t\t\t" << mzb->clusters[k].eta << endl;
		out_clusts << "   pT \t\t\t"  << mzb->clusters[k].pTtot << endl;
	    } //for each cluster
		
	    free(mzb->mcd);
	    free(mzb->clusters);
	    free(mzb);
	     
	} //while nevents < eventend
        out_clusts << "****" << nevents << " events****" << endl;
	free(tracks);
        in_tracks.close();
	out_clusts.close();
	plot1dat.close();
	plot2dat.close();
	ptdat.close();
	traxdat.close();
	disdat.close();
	cout << "5 new data files created." <<endl;
	return 0;
} //end main
