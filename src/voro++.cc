#include <voro++.hh>
#include <iostream>
#include <string>
#include <Gromacs.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <vector>
#include <fstream>
#include <thread>
#include <mutex>
#include <chrono>
#include <unistd.h>

using namespace std;
using namespace voro;

#define DEBUG

#define INPUT_TRJ_FILE 		1
#define INPUT_CONF_FILE 	2
#define OUTPUT_VORO_FILE 	3
#define COUNT_TREAD 		4
#define SET_BEGIN     		5
#define SET_END       		6
#define SET_STEP      		7

char* trj_file = "traj_comp.xtc";
char* conf_file = "conf.gro";
char* voro_file = "voro.vor";
mutex mtx;
md_file *mf; md_header mdh; md_ts mdts;
long int ts_first = 1; long int ts_step = 1; long int ts_last = LONG_MAX;
const int nbin =30;
int g_num_atom = 0;
int core = 1;
int ts_cnt=1; 
int comand_line( int argc, char** argv );
bool good_param( int i , int argc, char** argv);
int choise_flag( string choise );
bool check_read_header(md_file* mf , md_header* mdh);
string del_prob(string s);
bool read_norm();
void print_properties_step();
bool end_of_file_xdr();


class atom{
	public:
	int num_mol;
	string name_mol;
	string type_atom_in_mol;
	int num_atom = 0;
	float x=0,y=0,z=0;
	float vx=0,vy=0,vz=0;
	
	atom(string line);

};

atom::atom(string line){

	num_mol = stoi(line.substr(0,5));
	name_mol = del_prob(line.substr(5,5));
	type_atom_in_mol = del_prob(line.substr(10,5));
	g_num_atom++;
	num_atom=g_num_atom;
	x = stof(line.substr(20,8));
	y = stof(line.substr(28,8));
	z = stof(line.substr(36,8));

	if ( line.length() > 46 ){
		vx = stof(line.substr(44,8));
		vy = stof(line.substr(52,8));
		vz = stof(line.substr(60,8));
	};	

};

class molecul{
	public:
	int atom_to_mol=0;
};




class Point3{
	public:
	float	x;
	float	y;
	float  	z;
	float  	r;
};

class Frame{
	public:
		bool readOK;
		vector<Point3> pos;
		Point3 	       box;
	 
};

float get_r(char type){
	string type_atom;
	type_atom.push_back(type);
	const float rC = 0.175, rH = 0.125 , rO = 0.145 ;
	if(type_atom == "C" ){return rC;}
	if(type_atom == "O" ){return rO;}
	if(type_atom == "H" ){return rH;}
	return 0;
};



class grofile{
	public:
	int count_atom;
	int count_mol;
	vector<atom> all_atom;	
	grofile(string name_file);
	vector<float> r_list;

};

grofile::grofile(string name_file){
	string line;
	fstream f(name_file);
	if( !f.is_open() ){
		fprintf(stdout, "cannot open .gro file");
	} else {
		getline(f,line);
		getline(f,line);
		count_atom = stoi(line);
		for(int i = 0; i < count_atom; i++){
			getline(f,line);
			atom a(line);
			char atom_name = a.type_atom_in_mol[0];
			r_list.push_back(get_r(atom_name));
			all_atom.push_back(a);
		}
		f.close();
	}
};

class tess_voro{
	public: 
		vector<vector<int>> all_nei;
		vector<float> all_vol;

};


class Molecul{


	public:
	string molecul_type;
	vector<int> all_atom;
	vector<float> all_vol_atom;	
	float vol_mol;
	vector<int> nei_atom;
	vector<int> nei_mol;

}



Frame get_frame(grofile &gro);
void calc_voro_tess(const Frame &fr,int i, tess_voro &vor);

void build_molecul(tess_voro &vor, grofile gro){
	vector<Molecul> all_mol;
	for(int i = 0; i < vor.all_vol.size() ; i++){
		all_mol[gro.all_atom[i].num_mol].all_atom.push_back(i);
	}

	for(int i = 0;i < all_mol.size();i++){
                for(int j = 0; j < all_mol[i].all_atom.size(); j++){
			all_mol[i].all_vol_atom[j] = tess_voro.all_vol[all_mol[i].all_atom[j]]
		}
		float vol;
                for(int j = 0; j < all_mol[i].all_atom.size(); j++){
			vol+=all_mol[i].all_vol_atom[j];
		}

		all_mol[i].vol_mol = vol;

		for(int j = 0; j < all_mol[i].all_atom.size(); j++){
			for(int k = 0; k < vor.all_nei[all_mol[i].all_atom[j]].size(); k++){
				all_mol[i].nei_atom.push_back(vor.all_nei[all_mol[i].all_atom[j]][k]);
			}
		}

		for(int j = 0; j < ; j++){
		
        }

}


int main( int argc, char** argv ) {

	//КОМАНДНАЯ СТРОКА 
	comand_line( argc,argv );
	//OPEN_XTC_FILE
    	ts_last = 1;
    	/* Open Gromacs file with coordinates */
    	mf = mdio_open(trj_file, 0, MDIO_READ);
    	if(!mf) {fprintf(stderr, "%s: %s\n",mdio_errmsg(mdio_errno()), trj_file ); return 2;}
    	/* Read Gromacs header */
	if(!check_read_header(mf,&mdh)) return 3;
    	/* For .g96 files we don't have number of atoms in header */
    	if( mf->fmt == MDFMT_G96 ) mdh.natoms = g96_countatoms(mf);
    	mdts.natoms = mdh.natoms;
    	/* Read timesteps from file */
	
	grofile gro(conf_file);
	bool end = true;
	while(end){
		vector<Frame> all_fr;
		for( int i = 0; i < core; i++){
			Frame a;
		        a = get_frame(gro);
			if ( a.readOK == false ){ 
				end=false;break;
			} else {
			        all_fr.push_back(a);
			}
		}
		
		auto t1 = chrono::high_resolution_clock::now();

		vector<tess_voro> all_tess;	
		vector<thread> ths;
	
		for (int  i = 0; i < all_fr.size(); i++){
		       	all_tess.push_back(tess_voro());
		}

		for (int  i = 0; i < all_fr.size(); i++){
    			ths.push_back(thread( &calc_voro_tess, cref(all_fr[i]) , i , ref(all_tess[i])  ) );
		}
			
  		for (auto & th : ths)
    			th.join();


	    	auto t2 = chrono::high_resolution_clock::now();		
		auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    		chrono::duration<double, std::milli> ms_double = t2 - t1;
		
		cout << "time: "  << ms_int.count() << " ms\n";

		


/*		
		cout <<"num atom: "<< all_tess[0].all_vol.size() << endl;
		for(int i = 0 ; i < 10; i++){
			int size = all_tess[0].all_nei[i].size();
			cout <<"atom: "<< i <<" size:" << size <<" volume: " << all_tess[0].all_vol[i] <<"   nei: ";
			for(int j  = 0; j < size;j++){ 
				cout << all_tess[0].all_nei[i][j] << " ";}
			cout << endl;		

		}*/
		cout << "" <<endl;

    /* Getting number of milliseconds as a double. */
	}

}


string del_prob(string s){
	string v;
	for(char c:s) if (c != ' ') v += c;
	return v;
} 

int print_usage() {
    fprintf(stderr,"Usage: xtc2nice [-A] traj.[trr|xtc] [first last step]\n");
    exit(1);
}

int choise_flag( string choise ){
	if (choise == "-f"){return INPUT_TRJ_FILE   ;}
	if (choise == "-c"){return INPUT_CONF_FILE  ;}
	if (choise == "-o"){return OUTPUT_VORO_FILE ;}
	if (choise == "-j"){return COUNT_TREAD      ;}
	if (choise == "-b"){return SET_BEGIN        ;}
	if (choise == "-e"){return SET_END          ;}
	if (choise == "-s"){return SET_STEP         ;}

	return 0;

}
bool good_param( int i , int argc, char** argv){
	i++;
	if ( i < argc ){
		if ( choise_flag(argv[i]) ){ return false; }
		return true;
	}
	return false;
}

int comand_line( int argc, char** argv ){
	printf("argc: %d\n", argc);
	for(int i = 1; i < argc; i++){
		switch( choise_flag(argv[i]) ){
			case INPUT_TRJ_FILE: 
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;};
				i++;trj_file = argv[i];break;

			case INPUT_CONF_FILE: 
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;}; 
				i++;conf_file = argv[i];break;

			case OUTPUT_VORO_FILE: 
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;} 
				i++;voro_file = argv[i];break;

			case COUNT_TREAD: 
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;} 		
				i++;core = stoi(argv[i]);break;

			case SET_BEGIN:
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;} 		
				i++;ts_first = stoi(argv[i]);break;

			case SET_END  :
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;} 		
				i++;ts_last = stoi(argv[i]);break;

			case SET_STEP :
				if ( !good_param(i , argc , argv) ){cout << "error\n"; return 0;} 		
				i++;ts_step = stoi(argv[i]);break;

			case 0: cout << "not" << endl; break;
			
		}
	}

	cout << "comand line: -f "<< trj_file << " -c " << conf_file << " -o " << voro_file 
		<< " -j " << core << " -b " << ts_first << " -e " << ts_last << " -s " << ts_step << endl;
	return 0;
}



bool check_read_header(md_file* mf , md_header* mdh){
    	if( mdio_header(mf, mdh) < 0 ){
        	mdio_close(mf);
        	fprintf(stderr, "Cannot read header fromm '%s', %s\n",trj_file, mdio_errmsg(mdio_errno()));
        	return false;
    	}
	return true;
}



Frame get_frame(grofile &gro){
	Frame fr;
	if (!(mdio_timestep(mf, &mdts) < 0)){
		fr.box.x = mdts.box->A;
		fr.box.y = mdts.box->B;
		fr.box.z = mdts.box->C;
		fr.readOK = true;

		for(int i=0; i < mdts.natoms; ++i){
                	float *a = mdts.pos + 3*i;
			Point3 p;
			p.x = a[0];
			p.y = a[1];
			p.z = a[2];
			p.r = gro.r_list[i]; 
			fr.pos.push_back(p);
                }

		mdio_tsfree( &mdts );
		return fr;
	}else{
		fr.readOK = false;
		return fr;
	};

};

bool read_norm(){

	if( mdts.natoms != mdh.natoms ){
		fprintf(stderr, "Timestep in file contains wrong number of atoms!\n");
		fprintf(stderr, "Found %d, expected %d\n", mdts.natoms, mdh.natoms);
		return false;
	}

	if( !mdts.box ){
		fprintf(stderr, "Timestep does not contain box dimensions!\n");
		return false;
	}
	return true;
}

void print_properties_step(){
	fprintf(stdout, "#step: %-8d time: %g\n", mdts.step, mdts.time );
	fprintf(stdout,"#num atom: %d size box %g %g %g\n", 
		mdts.natoms , mdts.box->A ,mdts.box->B ,mdts.box->C );
}

bool end_of_file_xdr(){
	return (mdio_timestep(mf, &mdts) < 0);
}



void calc_voro_tess(const Frame &fr,int i,tess_voro &vor){
	container_poly con(0,fr.box.x ,0,fr.box.y,0,fr.box.z,nbin,nbin,nbin,true,true,true,4096);
	particle_order po;
	for(int  i=0; i < fr.pos.size(); i++){
		con.put(po , i , fr.pos[i].x , fr.pos[i].y , fr.pos[i].z , fr.pos[i].r);
	}
	voronoicell_neighbor c(con);
	c_loop_order vl(con, po);
	//int ind_mol = 0; 
	//tess_voro vor;


	auto t1 = chrono::high_resolution_clock::now();
	if(vl.start()){
		do{ 
			if(con.compute_cell(c,vl)) {

			//fprintf(stdout, "num atom %d; volume %e; neibors:", ind_mol , c.volume());
			vor.all_vol.push_back(c.volume());
			std::vector<int> v;
			c.neighbors(v);
			vor.all_nei.push_back(v);
			//for (auto it = v.begin(); it != v.end(); it++) {fprintf(stdout, "%d ", *it);}
			//fprintf(stdout, ";\n");
			}
			//ind_mol++;
		} while(vl.inc());
	}


	
    	auto t2 = chrono::high_resolution_clock::now();		
	auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
    	chrono::duration<double, std::milli> ms_double = t2 - t1;

	cout << i << " - " << this_thread::get_id() << "time: "  << ms_int.count() << "ms\n";
}
  /*


    	while( !end_of_file_xdr() )
    	{	
       		++ts_cnt;
		//if (read_norm()){return 0;};
		if( ts_cnt >= ts_first && ((ts_cnt - ts_first) % ts_step == 0) ){
			container con(0,mdts.box->A ,0,mdts.box->B,0,mdts.box->C,nbin,nbin,nbin,true,true,true,1024);
			particle_order po;
			for(int i=0; i < mdts.natoms; ++i){
				float *a = mdts.pos + 3*i;
				con.put(po,i,a[0],a[1],a[2]);
			}

			voronoicell_neighbor c(con);
			c_loop_order vl(con, po);
			int ind_mol = 0; 
			if(vl.start()){
			    	do{ 
					if(con.compute_cell(c,vl)) {
				   	fprintf(stdout, "num atom %d; volume %e; neibors:", ind_mol , c.volume());
				   	std::vector<int> v;
				   	c.neighbors(v);
				   	for (auto it = v.begin(); it != v.end(); it++) {
				      		fprintf(stdout, "%d ", *it);
				   	}
	
					fprintf(stdout, ";\n");
					}
					ind_mol++;
			    	} while(vl.inc());
			}
		}
		mdio_tsfree( &mdts );
		if( ts_cnt >= ts_last )break;
    	}

    	//fprintf(stderr, "ended read_timestep on step %d\n", ts_cnt);
//  */
