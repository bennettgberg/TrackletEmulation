#include "tracklet_em.h"

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
	int eventstart = 1;
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
	int nevents = 0;
	getline(in_tracks, data_in, ' ');
        string data;
 //pointer to structure to hold jet MC data for each event.
        struct mc_data * mcdat = (struct mc_data *)malloc(sizeof(struct mc_data));
     	while(nevents < eventend){
		int ntracks = 0;
		if(data_in == "Event"){
			++nevents;
			for(int cc=0; cc<4; ++cc){
				getline(in_tracks, data_in, ' ');
			}
			getline(in_tracks, data, ' ');
                        mcdat->ogpt = atof(data.c_str());
			for(int cc=0; cc<3; ++cc){
				getline(in_tracks, data, ' ');
			}
                        mcdat->ogeta = atof(data.c_str());
			getline(in_tracks, data);
                        data = data.substr(0, data.length()-1);
                        mcdat->ogphi = atof(data.c_str());
		      //Next line gives no important information.
			getline(in_tracks, data_in);
			if(nevents < eventstart){
				ntracks = 0;
				continue;
			}
	   		out_clusts << "\n****EVENT " << nevents << " ****" << endl;
	   		out_clusts << "pT: " << mcdat->ogpt << " eta: " << mcdat->ogeta << " phi: " << mcdat->ogphi << endl;
			getline(in_tracks, data_in, ' ');
		}
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
             mcdat->ntracks = ntracks;
             struct maxzbin * mzb = L2_cluster(tracks, mcdat, nzbins);
             if(mzb == NULL){
                continue;
             }
         //output to plot data files.
             for(int k = 0; k < mzb->nclust; ++k){
               //if jet and cluster phi and eta are opposite sign, cluster belongs to opposite jet.
		if (mzb->phimc*mzb->clusters[k].phi<0){
			if(mzb->phimc<0)
				plot1dat << mzb->phimc+M_PI << "\t" << mzb->clusters[k].phi << endl;
			else
				plot1dat << mzb->phimc-M_PI << "\t" << mzb->clusters[k].phi << endl;
		}
		else
	     		plot1dat << mzb->phimc << "\t" << mzb->clusters[k].phi << endl;
		if (mzb->etamc*mzb->clusters[k].eta<0)
			plot2dat << mzb->etamc*-1 << "\t" << mzb->clusters[k].eta << endl;
		else
	        	plot2dat << mzb->etamc << "\t" << mzb->clusters[k].eta << endl;
		ptdat << mzb->pTmc << "\t" << mzb->clusters[k].pTtot << endl;
		traxdat << ntracks << "\t" <<  mzb->clusters[k].numtracks << endl;
	     }
		
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
	cout << "4 new data files created." <<endl;
	return 0;
} //end main
