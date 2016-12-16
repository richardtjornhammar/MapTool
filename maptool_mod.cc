//// C/C++ STUFF
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>

//// CLIPPER STUFF
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

//// GSL STUFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define	DIM	3

// CHECK: http://www.ccp4.ac.uk/html/edstats.html
//-----------------------------------------
//	OLD COMPILE USING:
// export LD_LIBRARY_PATH=/usr/local/lib && g++ main.cc -lclipper-core -lclipper-contrib -lclipper-ccp4 -lclipper-phs -o richmap
//-----------------------------------------
//	COMPILE:
// export LD_LIBRARY_PATH=/home/richard/autobuild/Linux-Carmack-pre-release-gtk2-noguile/lib && g++ maptool.cc -lclipper-core -lclipper-contrib -lclipper-ccp4 -lclipper-phs -lgsl -lblas -o richmap
//-----------------------------------------

class MapFilterFn_g5 : public clipper::MapFilterFn_base {
	public: 
		clipper::ftype operator() ( const clipper::ftype& radius ) const {
			return exp(-radius*radius/50.0); 
		}
};


void
fatal(void){
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" 	<< std::endl;
	std::cout << "PLEASE CONTACT "		<< author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

void
fatal( std::string errstr ) {
	std::string author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" 	<< std::endl;
	std::cout << "INFO:: " 		<< errstr << std::endl;
	std::cout << "PLEASE CONTACT " 	<< author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

void
die( std::string errstr ){
	fatal( errstr );
}


std::pair< int , std::vector< std::string > >
argparser(	std::pair	< int, char ** > mp, std::vector<int>& v_set,
		std::vector	< std::pair < int, std::pair<std::string, std::string > > > opts_n_defaults ) {
	//
	// ARGS OUTPUT ORDER CORRESPONDS TO OPTS_N_DEFAULTS FIRST VALUE
	//
	std::vector< std::string > args;
	std::pair  < int, std::vector< std::string > > ret_args;
	//
	//	std::vector< int > v_set;
	//
	for( int i=0 ; i<mp.first ; i++ ) {
		args.push_back			( mp.second[i] );
	}
	for( int i=0 ; i<opts_n_defaults.size() ; i++ ) {
		ret_args.second.push_back	( opts_n_defaults[i].second.second );
		// v_set.push_back(0);
	}
	// COMMAND INPUT
	int arg 	= 0 ;
	ret_args.first 	= 0 ;
	int a_size = mp.first==1?2:mp.first;
	int help=0,failed=0;
	while ( ++arg < a_size  ) {
		if(  args[arg] == "-h" ||  args[arg] =="--h" || args[arg] =="--help" || args[arg] =="-help" || mp.first==1) {
			help=1;
			break;
		}else{
			int ia=arg,iap=++arg;
			for(int i=0;i<opts_n_defaults.size();i++) {
				if ( args[ia] == opts_n_defaults[i].second.first ) {
					ret_args.second[i] = args[iap];
					v_set[i]=1;
				}
			}
		}
	}
	for(int i=0;i<v_set.size();i++) {
		if( opts_n_defaults[i].first > 0 && v_set[i] == 0 ) {
			ret_args.first = 1;
			failed=1;
			break;
		}
	}
	if ( help==1 || failed==1 ) {
		std::cout << "HELP:: VIABLE INPUT OPTIONS ARE:\nHELP:: [X] \t OPTION \t VALUE\n" << std::endl;
		for(int i=0;i<opts_n_defaults.size();i++)
			std::cout << "HELP:: [" << opts_n_defaults[i].first << "] \t "
				<< opts_n_defaults[i].second.first << " \t "
				<< ret_args.second[i] << std::endl;
		ret_args.first = 1;
		std::cout << std::endl;
		if(failed==1) {
			std::cout << "ERROR::OBLIGATORY ARGUMENTS NOT SET ( MARKED [1] )\n" << std::endl;
		}
	}

	return ret_args;
}


void SplitFilename (const std::string& str)
{
  std::cout << "Splitting: " << str << '\n';
  std::size_t found = str.find_last_of("\\/");
  std::cout << " path: " << str.substr(0,found) << '\n';
  std::cout << " file: " << str.substr(found+1) << '\n';
}

std::string SplitString (const std::string& str,const std::string& split_str, int i)
{
	std::size_t found = str.find_last_of(split_str);
	std::string retstr;
	switch( i ) {
		case 1:
			retstr = str.substr(found+1);
			break;
		default:
			retstr = str.substr(0,found);
	}
	return retstr;	
}

clipper::Xmap<float> apply_filter(clipper::Xmap<float> base_map, float filter[3][3][3]) {

	clipper::Xmap_base::Map_reference_index midx = base_map.first();
	clipper::Xmap<float> f_map(base_map);
	float cnt = 27.0;
	int I = 1;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		float vox_ave=0.0;
		for(int i=-1;i<=1;i++){
			for(int j=-1;j<=1;j++){
				for(int k=-1;k<=1;k++){
					vox_ave	+= filter[i+1][j+1][k+1]*base_map.get_data( midx.index_offset(i,j,k) );
				}
			}
		}
		float val = vox_ave/cnt;
		f_map.set_data( midx.index(), val );
		I++;
	}

	return f_map;
}

clipper::Xmap<float> unsharp_masking(clipper::Xmap<float> base_map ){

	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	clipper::Xmap<float> um_map(base_map);
	int I = 1;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		float vox_ave=0.0, cnt=0.0;
		for(int i=-1;i<=1;i++){
			for(int j=-1;j<=1;j++){
				for(int k=-1;k<=1;k++){
					vox_ave	+= base_map.get_data( midx.index_offset(i,j,k) );
					cnt	+= 1.0;	
					// std::cout << "INFO:: " << i << ", " << j << ", " << ", " << k << std::endl;
				}
			}
		}
		float val = base_map.get_data( midx.index() ) - vox_ave/cnt;
		val	  = val>0?val:0;
		um_map.set_data( midx.index(), val );
		I++;
	}
	return um_map;
}

clipper::Xmap<float> blur_map(clipper::Xmap<float> base_map ){ // index offset only works for +-1
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	clipper::Xmap<float> b_map(base_map);
	int I = 1;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		float vox_ave=0.0, cnt=0.0;
		for(int i=-1;i<=1;i++){
			for(int j=-1;j<=1;j++){
				for(int k=-1;k<=1;k++){
					vox_ave	+= base_map.get_data( midx.index_offset(i,j,k) );
					cnt	+= 1.0;	
				}
			}
		}
		float val = vox_ave/cnt;
		val	  = val>0?val:0;
		b_map.set_data( midx.index(), val );
		I++;
	}
	return b_map;
}

int sgn(int i){
	return i>=0?1:-1;
}

clipper::Xmap<float> sobel_map(clipper::Xmap<float> base_map ){ 
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	clipper::Xmap<float> b_map(base_map);
	int I = 1;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		float vox_ave=0.0, cnt=0.0;
		// ONE FOR EACH  u,v,w
		float d[3]={1,2,1};

		for(int i=-1;i<=1;i++){
			for(int j=-1;j<=1;j++){
				for(int k=-1;k<=1;k++){
					vox_ave	+= base_map.get_data( midx.index_offset(i,j,k))*d[i+1]*j;
					vox_ave	+= base_map.get_data( midx.index_offset(i,j,k))*d[j+1]*k;
					vox_ave	+= base_map.get_data( midx.index_offset(i,j,k))*d[k+1]*i;
					cnt	+= 3.0;	
				}
			}
		}

		float val = vox_ave/cnt;
		val	  = val>0?val:0;
		b_map.set_data( midx.index(), val );
		I++;
	}
	return b_map;
}

clipper::Xmap<float> random_blur_map(clipper::Xmap<float> base_map ) { 
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	clipper::Xmap<float> b_map(base_map);

//	RANDOM NUMBER GENERATION
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		float vox_ave = 0.0, cnt = 0.0;
		for(int i=-1;i<=1;i++) {
			for(int j=-1;j<=1;j++) {
				for(int k=-1;k<=1;k++) {
					float r_val	 = gsl_ran_flat ( r, 0.0, 1.0); // average will be 1/2
					vox_ave		+= base_map.get_data( midx.index_offset(i,j,k))*r_val;
					cnt		+= r_val;
				}
			}
		}
		float val = 2.0*vox_ave/cnt;
		b_map.set_data( midx.index(), val );
	}

	gsl_rng_free(r);

	return b_map;
}

float calc_kurtosis(clipper::Xmap<float> base_map, int verbose) {
//
//	DOING MAP STATISTICS, RETURNS KURTOSIS	
//
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	float calc_fs[5]={0.0, 0.0, 0.0, 0.0, 0.0}, e_m, e_s2;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		calc_fs[0] += base_map[midx] ;
		calc_fs[1] += base_map[midx]*base_map[midx] ;
		calc_fs[2] += 1.0;
	}
	e_m 		= calc_fs[0]/calc_fs[2];
	calc_fs[3]	= calc_fs[1]/calc_fs[2];
	e_s2		= calc_fs[3]-e_m;
	if(verbose)
		std::cout << "STAT::" << e_m << std::endl;
	double cookie = 0.0, nom = 0.0, nam = 0.0;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		nom += (base_map[midx]-e_m)*(base_map[midx]-e_m)*(base_map[midx]-e_m)*(base_map[midx]-e_m);
		nam += (base_map[midx]-e_m)*(base_map[midx]-e_m);
	}

	if(verbose>1)
		std::cout << "STAT:: " << nom << ", " << nam << ", " << calc_fs[2] << std::endl;
	cookie   = nom/nam/nam*calc_fs[2];
	return cookie;
}

float calc_skewness(clipper::Xmap<float> base_map, int verbose) {
//
//	DOING MAP STATISTICS, RETURNS SKEWNESS	
//
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	float calc_fs[5]={0.0, 0.0, 0.0, 0.0, 0.0}, e_m, e_s2;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		calc_fs[0] += base_map[midx] ;
		calc_fs[1] += base_map[midx]*base_map[midx] ;
		calc_fs[2] += 1.0;
	}
	e_m 		= calc_fs[0]/calc_fs[2];
	calc_fs[3]	= calc_fs[1]/calc_fs[2];
	e_s2		= calc_fs[3]-e_m;
	if(verbose)
		std::cout << "STAT::" << e_m << std::endl;
	double cookie = 0.0, nom = 0.0, nam = 0.0;
	for(midx = base_map.first(); !midx.last(); midx.next() ) {
		nom += (base_map[midx]-e_m)*(base_map[midx]-e_m)*(base_map[midx]-e_m);
		nam += (base_map[midx]-e_m)*(base_map[midx]-e_m);
	}

	if(verbose>1)
		std::cout << "STAT:: " << nom << ", " << nam << ", " << calc_fs[2] << std::endl;
	double	om	= sqrt(nam);
	double	namnam	= calc_fs[2];
	
	cookie   = nom/nam/om*namnam;
	return cookie;
}


float calc_standardized_moment(clipper::Xmap<float> base_map, int imom , int verbose) {
//
//	DOING MAP STATISTICS, RETURNS STANDARDIZED MOMENT
//
	clipper::Xmap_base::Map_reference_index midx=base_map.first();
	float calc_fs[5]={0.0, 0.0, 0.0, 0.0, 0.0}, e_m, e_s2;
	for( midx = base_map.first(); !midx.last(); midx.next() ) {
		calc_fs[0] += base_map[midx] ;
		calc_fs[1] += base_map[midx]*base_map[midx] ;
		calc_fs[2] += 1.0;
	}
	e_m 		= calc_fs[0]/calc_fs[2];
	calc_fs[3]	= calc_fs[1]/calc_fs[2];
	e_s2		= calc_fs[3]-e_m;
	if(verbose)
		std::cout << "STAT::" << e_m << std::endl;
	double cookie = 0.0, nom = 0.0, nam = 0.0;
	for( midx = base_map.first() ; !midx.last() ; midx.next() ) {
		double nmomv = (base_map[midx]-e_m)*(base_map[midx]-e_m);
		for( int im=0 ; im<imom ; im++ )
			nmomv*=(base_map[midx]-e_m);
		nom +=  nmomv;
		nam += (base_map[midx]-e_m)*(base_map[midx]-e_m);
	}
	if(verbose>1)
		std::cout << "STAT:: " << nom << ", " << nam << ", " << calc_fs[2] << std::endl;
	double	namnam	= calc_fs[2];	
	nam = pow(nam,imom+2.0);
	double	om	= sqrt(nam);
	cookie   	= nom/om*namnam;
	return cookie;
}


float calc_map_mean( clipper::Xmap<float> base_map , int bAbs) {
	float nom=0.0, nam=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = base_map.first(); !midx.last(); midx.next() ) {
		nom += bAbs>0?abs(base_map[midx]):base_map[midx];
		nam += 1.0;
	}
	float e_m=nom/nam;
	return e_m;
}

float sum_map( clipper::Xmap<float> base_map , int bAbs) {
	float sum=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = base_map.first(); !midx.last(); midx.next() ) {
		sum += bAbs>0?abs(base_map[midx]):base_map[midx];
	}
	return sum;
}

clipper::Xmap<float> calc_map_scale( clipper::Xmap<float> map1 , float scale) {
	clipper::Xmap<float> map3(map1);
	for(clipper::Xmap<float>::Map_reference_index midx = map1.first(); !midx.last(); midx.next() ) {
		map3[midx]=map1[midx]*scale;
	}
	return map3;
}

clipper::Xmap<float> calc_map_mult( clipper::Xmap<float> map1 , clipper::Xmap<float> map2) {
	float nom=0.0, nam=0.0;
	clipper::Xmap<float> map3(map1);
	for(clipper::Xmap<float>::Map_reference_index midx = map2.first(); !midx.last(); midx.next() ) {
		map3[midx]=map1[midx]*map2[midx];
	}
	return map3;
}

std::string conv_f2s (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();   
}

std::vector<float > calc_map_extreme( clipper::Xmap<float> base_map) {
	float maximum=0.0, minimum=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = base_map.first(); !midx.last(); midx.next() ) {
		maximum=(base_map[midx]>maximum)?(base_map[midx]):(maximum);
		minimum=(base_map[midx]<minimum)?(base_map[midx]):(minimum);
	}
	std::vector<float > extremes;
	extremes.push_back(maximum);	
	extremes.push_back(minimum);
	return extremes;
}

clipper::Xmap<float> histogram_equalisation( clipper::Xmap<float> map , std::vector<std::pair<float, float[2] > > stats) {
	clipper::Xmap<float> map_out(map);
	float lims[2]={stats[0].second[0],stats[0].second[1]}, maximal = lims[0]-lims[1];
	int   Nb = ((int)stats[0].first);
	std::cout << "INFO::HE:: "<< lims[0] << " " << lims[1] << " " << Nb << " " << maximal << std::endl;

	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		float val	= map[midx]-lims[1];
		int I		= ceil( val*stats[0].first/maximal );  
		I		= ( (I<1)?(1):((I>Nb)?(Nb):(I)) ); 
		float heval	= stats[I].second[1]*maximal+lims[1];
		map_out[midx]	= heval;
	}

	return map_out;
}

clipper::Xmap<float> stretch_equalisation( clipper::Xmap<float> map ,
					   std::vector<std::pair<float, float[2] > > stats,
					   float cutoff1, float cutoff2 ) 
{
	clipper::Xmap<float> map_out(map);
	float lims[2]={stats[0].second[0],stats[0].second[1]}, maximal = lims[0]-lims[1];
	int   Nb = ((int)stats[0].first);
	float d=lims[0], c=lims[1], a=0.0, b=0.0;
	for( int i = 1; i <= Nb ; i++ ){
		a = (stats[i].second[1]<cutoff1)?(stats[i].second[1]):(a);
		b = (stats[i].second[1]<cutoff2)?(stats[i].second[1]):(b);
	}
	std::cout << "INFO:: HAVE a,b,c,d = " << a << " " << b << " " << c << " " << d << std::endl;

	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		float pout 	= ( map[midx] - c )*( ( b - a )/( d - c ) ) + a ;
		map_out[midx]	= ( pout>a && pout<b ) ? ( pout ):( 0.0 );
	}

	return map_out;
}

clipper::Xmap<float> equalisation( clipper::Xmap<float> map ,
				std::vector<std::pair<float, float[2] > > stats,
				int type ) {
	clipper::Xmap<float> map_out(map);
	float lims[2]={stats[0].second[0],stats[0].second[1]}, maximal = lims[0]-lims[1];
	int   Nb = ((int)stats[0].first);

	if( type != 0 && type != 1 ){
		std::cout << "INFO::ERROR:: "<< lims[0] << " " << lims[1] << " " << Nb << " " << maximal << std::endl;
		return map;
	}

	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		float val	= map[midx]-lims[1];
		int I		= ceil( val*stats[0].first/maximal ); // NOTE: 1->Nb
		I		= ( (I<1)?(1):((I>Nb)?(Nb):(I)) );
		float heval	= stats[I].second[type]*val+lims[1]; //stats[I].second[type]*maximal+lims[1];
		map_out[midx]	= heval;
	}

	return map_out;
}

void output_stats(std::vector<std::pair<float, float[2] > > stats, std::string fname) {
	float lims[2]={stats[0].second[0],stats[0].second[1]}, maximal = lims[0]-lims[1], di=maximal/stats[0].first;
	int   Nb = ((int)stats[0].first);
  	std::ofstream ofs;
  	ofs.open(fname.c_str()); //, std::ofstream::out | std::ofstream::app);
	for( int i = 1; i <= Nb; i++ ) {
		ofs << i*di+lims[1] << "\t" << stats[i].second[0] << "\t" << stats[i].second[1] << std::endl;
	}
	ofs.close();
}

std::vector<std::pair<float, float[2] > >
calc_ideal_peak_cdf( std::vector<std::pair<float, float[2] > > stats, std::string fname ) {
	float lims[2]	= {stats[0].second[0],stats[0].second[1]}, maximal = lims[0]-lims[1], di=maximal/stats[0].first;
	int   Nb	= ((int)stats[0].first);
	float dx	= (lims[0] - lims[1])/stats[0].first;
	float sigma	= stats[Nb+1].second[0], mu = stats[Nb+1].first;
	float norm	= gsl_sf_erf((stats[Nb].first-mu)/sigma)*2.0, zero = gsl_sf_erf((stats[1].first-mu)/sigma);
	float diff_val  = 0.0, diff_norm=0.0;
	for( int i = 1; i <= Nb ; i++ ) {
		float val  = ( gsl_sf_erf((stats[i].first-mu)/sigma)-zero )/norm; 	// AN IDEAL CDF
		diff_val   = val-stats[i].second[1];					// THE DERIVATIVE (PROBABILITY)
		stats[i].second[0] = diff_val;
		diff_norm += stats[i].second[0];
	}

////	NORMALIZE PRUNED PROBABILITY AND CDF
	for( int i = 1; i <= Nb ; i++ ) {
		stats[i].second[0] /= diff_norm*dx;
		stats[i].second[1]  = 0.0;
		for( int j=i; j>=1 ; j-- )
			stats[i].second[1] += stats[j].second[0]*dx;
	}

////	FILE IO
  	std::ofstream ofs;
  	ofs.open(fname.c_str());
	for( int i = 1; i <= Nb ; i++ )
		ofs << stats[i].first << "\t" << stats[i].second[0] << "\t" << stats[i].second[1] << std::endl;
	ofs.close();

	return stats;
}

int convolve_xmap( clipper::Xmap<float> manip_map )
{
	clipper::Skeleton_basic::Neighbours neighb(manip_map);
	int 	n_neighbs	= neighb.size();
	float	f_neig		= neighb.size();	
	float	v, v0;

	clipper::Xmap_base::Map_reference_index	ix;
	clipper::Coord_grid 			c_g;
	clipper::Coord_grid			c_g_start;

//	BLURRING
	for (ix = manip_map.first(); !ix.last(); ix.next()) {
		v0		 = manip_map[ix];
		c_g_start	 = ix.coord();
		for (int i=0; i<n_neighbs; i++) {
			c_g	 = c_g_start + neighb[i];
			v	+= manip_map.get_data(c_g);
		}
		manip_map[ix] += v/f_neig; 
	}

	return 0;
}

std::pair<float,float>
map_mean_variance(clipper::Xmap<float> map, bool bAbs){
	std::pair< float, float > mv;
	float X=0.0,X2=0.0,C=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		float val 	= map[midx];
		val		= (bAbs)?(sqrt(val*val)):(val);
		X		+=val;
		X2		+=val*val;
		C+=1.0;
	}
	X/=C; X2/=C;
	mv.first=X; mv.second=sqrt(X2-X*X);
	return mv;
}

std::pair< float, float  >
calc_map_mean_std( clipper::Xmap<float> map ) 
{
	float mean=0.0	, stdev=0.0, X=0.0, X2=0.0, C=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		X		+=map[midx];
		X2		+=map[midx]*map[midx];
		C++;
	}
	mean=X/C; X2/=C; 
	stdev=sqrt(X2-X*X);
	std::pair< float, float  > ms;
	ms.first=mean; ms.second=stdev;
	return ms;
}

std::vector<std::pair<float, float[2] > >
calc_cdf(clipper::Xmap<float> map, std::vector<float > lims) {
	std::vector<float > cdf, prob;
	float range = lims[0] - lims[1];
	int Nb = 512;
//FOR ENTIRE MAPSTATS
	float mean=0.0, stdev=0.0, X=0.0, X2=0.0, C=0.0;
//B FOR FWHM CALC
	float pmax=0.0, hmax=0.0, sh[3]={0.0,0.0,0.0};
	int   imax=0, ihmax=0;
//E FFC
	for(int i=0;i<Nb;i++){
		cdf .push_back(0.0);
		prob.push_back(0.0);
	}
	float prob_sum=0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = map.first(); !midx.last(); midx.next() ) {
		float val 	= map[midx] - lims[1];
		X		+=map[midx];
		X2		+=map[midx]*map[midx];
		C++;
		int I		= floor( ((float)Nb)*val/range );
		if( I>=0 && I<Nb ){
			prob[I]		+= 1.0;
			prob_sum	+= 1.0;
		}
	}
//	OF ENTIRE DISTRIBUTION
	mean=X/C; X2/=C; stdev=sqrt(X2-X*X);
	std::cout << "INFO::STAT ::M:: " << mean << " ::S:: " << stdev << " | " << C << std::endl;
	std::pair<float, float[2] > stat_pack;
	std::vector<std::pair<float, float[2] > > stats;
	stat_pack.first	= ((float)Nb);
	stat_pack.second[0]	= lims[0];
	stat_pack.second[1]	= lims[1];
	// std::cout << "INFO::LIMITS " << lims[0] << " :::: " << lims[1] << std::endl;
	stats.push_back(stat_pack);
	for ( int i=0 ; i<Nb ; i++ ) {
		prob[i]/=prob_sum;
		for( int j=i; j>=0 ; j-- )
			cdf[i] += prob[j];
		stat_pack.first		= i*range/((float)Nb)+lims[1];
		stat_pack.second[0]	= prob[i]*((float)Nb)/(lims[0]-lims[1]);
		stat_pack.second[1]	= cdf[i];
//B FOR FWHM CALC
		imax  = (stat_pack.second[0]>=pmax)?(i):(imax);
		sh[0] = (imax==i)?(stat_pack.first):(sh[0]);
		pmax  = (stat_pack.second[0]>=pmax)?(stat_pack.second[0]):(pmax);
		hmax  = 0.5*pmax;
		ihmax = (stat_pack.second[0]>=0.5*pmax)?(i):(ihmax);
		sh[1] = (ihmax==i)?(stat_pack.first):(sh[1]);
		hmax  = (stat_pack.second[0]>=0.5*pmax)?(stat_pack.second[0]):(hmax);
//E FFC
		stats.push_back(stat_pack);
  	}
	for ( int i=1 ; i<=Nb ; i++ ) {
		if( i<imax ) {
			sh[1]=(stats[i].second[0]<hmax)?stats[i].first:sh[1];
		} else {
			sh[2]=(stats[i].second[0]>hmax)?stats[i].first:sh[2];
		}
	}

	//// THE PEAK PROPERTIES ////
	// std::cout << "INFO::STATP::PM:: " << pmax << " ::HM:: " << hmax << std::endl;	
	// std::cout << "INFO::STATP::PP:: " << stats[imax].first << std::endl;
	// std::cout << "INFO::STATP::M::  " << (sh[2]+sh[1])*0.5 << " ::S:: " << (sh[2]-sh[1])*0.5 << std::endl;	
	stat_pack.first		= stats[imax].first; // (sh[2]+sh[1]) * 0.5;
	stat_pack.second[0]	= (sh[2]-sh[1])*0.5;
	stat_pack.second[1]	=  sh[2];
	stats.push_back(stat_pack);

	return stats;
}

// REMEMBER
// g(x)=1/sqrt(2*pi*sigma^2)*exp( -0.5*( (x-mu)^2/sigma^2 ) )
//
// integral becomes
// sqrt(\sigma)*erf( \frac{ x - \mu }{ sqrt(2) \sigma } ) + 2 * \sigma


clipper::Xmap<float> quench_map( clipper::Xmap<float> map0 ) {
	clipper::Xmap<float> map3(map0);
	float x = 0.0;
	for(clipper::Xmap<float>::Map_reference_index midx = map0.first(); !midx.last(); midx.next() ) {
		x=map0[midx];
		map3[midx]=x*(tanh((x-0.4)*4.0*M_PI)+1.0)/2.0; // tanh is slightly faster than erf
	}
	return map3;
}

clipper::Xmap<float> degaus( clipper::Xmap<float> I_in, float s1, float s2 , bool full){

	clipper::Xmap<float> I_out(I_in);

	float s1_2 	= s1*s1; 
	float s2_2 	= s2*s2;
	float x 	= 0.0;

	if(full){ // THIS METHOD IS REALLY BAD SINCE IT LOOPS OVER EVERYTHING
		for( clipper::Xmap<float>::Map_reference_index iidx = I_in.first(); !iidx.last(); iidx.next() ) {
			float val		= 0.0;
			float cnt		= 0.0;
			float Ir		= I_in[iidx];
			clipper::Coord_grid c_r	= iidx.coord();
			for( clipper::Xmap<float>::Map_reference_index jidx = I_in.first(); !jidx.last(); jidx.next() ) {
				float Ix		= I_in[jidx];
				clipper::Coord_grid c_x	= jidx.coord();
				clipper::Coord_grid c_d; 
				c_d 		 = c_r - c_x;
				float du2 	 = c_d.u()*c_d.u() + c_d.v()*c_d.v() + c_d.w()*c_d.w();
				float dI2 	 = (Ir-Ix)*(Ir-Ix);
				float w		 = exp(-0.5*du2/s1_2)*exp(-0.5*dI2/s2_2);
				val		+= w*Ix;
				cnt		+= w;
			}
			I_out[iidx]=val/cnt;
		}
	}else{
		clipper::Skeleton_basic::Neighbours neighb(I_in, 0.1, s1*4.0 );
		int 	n_neighbs	= neighb.size();
		float	f_neig		= neighb.size();	
		float	Ir, Ix; 				
		float   val		= 0.0;
		clipper::Xmap_base::Map_reference_index	ix;
		clipper::Coord_grid 			c_r,c_x,c_d;

		for ( ix = I_in.first(); !ix.last(); ix.next() ) {
			Ir	 = I_in[ix];
			c_r	 = ix.coord();
			val	 = 0.0;
			float ws	= 0.0;
			for (int i = 0 ; i<n_neighbs ; i++ ) {
				c_x	 	 = c_r + neighb[i];
				c_d	 	 = neighb[i];
				Ix  	 	 = I_in.get_data(c_x);
				float du2 	 = c_d.u()*c_d.u() + c_d.v()*c_d.v() + c_d.w()*c_d.w();
				float dI2 	 = (Ir-Ix)*(Ir-Ix);
				float w		 = exp(-0.5*du2/s1_2)*exp(-0.5*dI2/s2_2);
				val		+= w*Ix;
				ws		+= w;
				
			}
			I_out[ix] = val/ws; 
		}
	}

	return I_out;
}


clipper::Xmap<float> 
degaus( clipper::Xmap<float> I_in, std::vector< float > stg ) 
{
	clipper::Xmap<float> I_out( I_in );

	float x 	= 0.0 ;
	float kurt	= calc_kurtosis( I_in, 0 );	
	std::pair< float , float > ms = map_mean_variance(I_in,false);

	std::cout << "::INFO::KURT::INFO::" << kurt << std::endl;
	std::cout << "::INFO::STAT::INFO::" << ms.first	 << std::endl;
	std::cout << "::INFO::STAT::INFO::" << ms.second << std::endl;

	float n2pi = sqrt(2.0*M_PI);

	float s1_2 = (stg.size()==3) ? (stg[0]*stg[0]) : ( kurt );
	float s2_2 = (stg.size()==3) ? (stg[1]*stg[1]) : ( ms.second*ms.second*2.0 );
	float s3_2 = (stg.size()==3) ? (stg[2]*stg[2]) : ( ms.second*ms.second*4.0 );

	float iG_1 = 1.0/( n2pi*sqrt(s1_2) );
	float iG_2 = 1.0/( n2pi*sqrt(s2_2) );
	float iG_3 = 1.0/( n2pi*sqrt(s3_2) );
/*
	kurt		= 1.0*kurt
	sigmaI1	 	= 1.0*np.std(imap)**2
	sigmaI2	 	= 1.0*np.std(imap)**2
	sigmap		= np.sqrt(kurt)**2
	npix 		= int( round(kurt*5) )

	n_dr_map	= imap_o
	sq2g	= np.sqrt(2.0*np.pi)
	iG1 	= 1.0/(sq2g*np.sqrt(sigmaI1))
	iG2 	= 1.0/(sq2g*np.sqrt(sigmaI2))
	iG3 	= 1.0/(sq2g*np.sqrt(sigmap ))
*/
	float s1 = sqrt(s1_2);
	clipper::Skeleton_basic::Neighbours neighb( I_in, 0.1, s1*5.0 );
	int 	n_neighbs	= neighb.size();
	float	f_neig		= neighb.size();	
	float	Ir, Ix; 				
	float   val		= 0.0;
	clipper::Xmap_base::Map_reference_index	ix;
	clipper::Coord_grid 	c_r,c_x,c_d;
	for ( ix = I_in.first(); !ix.last(); ix.next() ) {
		Ir		= I_in[ix];
		c_r		= ix.coord();
		val	 	= 0.0;
		float ws	= 0.0;
		for (int i = 0 ; i<n_neighbs ; i++ ) {
			c_x	 	 = c_r + neighb[i];
			c_d	 	 = neighb[i];
			Ix  	 	 = I_in.get_data(c_x);
			float du2 	 = c_d.u()*c_d.u() + c_d.v()*c_d.v() + c_d.w()*c_d.w();
			float dI2 	 = (Ir-Ix)*(Ir-Ix);
			float I2 	 = (Ir*Ix);
			float w		 = exp(-0.5*du2/s1_2)*iG_1*exp(-0.5*dI2/s2_2)*iG_2*exp(-0.5*I2/s3_2)*iG_3;
			val		+= w*Ix;
			ws		+= w;			
		}
		I_out[ix] = val/ws; 
	}

	return I_out;
}

clipper::Xmap<float> quench_map_with_cdf( clipper::Xmap<float> map0 ) {
	clipper::Xmap<float> map3( map0 );

	std::vector<float > lims	= calc_map_extreme( map3 );
	std::vector<std::pair<float, float[2] > > stats;
	stats				= calc_cdf( map3 , lims );
	float maximal	= lims[0]-lims[1];
	int   Nb 	= ((int)stats[0].first);
	float ppos	= stats[ stats.size()-1 ].first;
	float spread	= stats[ stats.size()-1 ].second[0];
	float hpos	= stats[ stats.size()-1 ].second[1];
	for (int i=0;i<stats.size()-1;i++) 	// use cdf 0.9% position
		if( stats[i].second[1] < 0.9 )
			hpos=stats[i].first;
	//std::cout << "INFO::QUENCH::STATS:: " << ppos << " " << hpos << " " << spread << std::endl;
	float x 	= 0.0;
	float fmap0	= sum_map( map0 , 1 );
	for( clipper::Xmap<float>::Map_reference_index midx = map0.first(); !midx.last(); midx.next() ) {
		x		 = map0[midx];					// THE VALUE
		float val	 = x-lims[1];		
		int I		 = ceil( val*stats[0].first/maximal ); 
		I		 = ( (I<1)?(1):((I>Nb)?(Nb):(I)) ); 
		float cdf_val	 = stats[I].second[1];				// THE CDF
		map3[midx]	 = x;
		float s 	 = ( erf((x-hpos)/spread)+1.0 )/2.0; 
		map3[midx]	*= s;
	}
	float fmap3	= sum_map( map3 , 1 );
	float tot_scale	= fmap0/fmap3;
	//std::cout << "INFO::SCALING::" << tot_scale  << std::endl;
	for( clipper::Xmap<float>::Map_reference_index midx = map0.first(); !midx.last(); midx.next() ) {
		map3[midx]	*= tot_scale;
	}
	return map3;
}


clipper::Xmap<float>
sharpen ( float b_factor, clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis ) {

	int n_data = 0;
	int n_inf  = 0;

	clipper::Grid_sampling gs(	original_fphis.spacegroup(), 
					original_fphis.cell(), 
					original_fphis.resolution(), 1.5 );
	clipper::Xmap<float> xmap;	  					// NULL MAP
	xmap.init ( original_fphis.spacegroup(), original_fphis.cell(), gs ); 	// INITIALIZE MAP

	clipper::HKL_info::HKL_reference_index hri;
	clipper::HKL_data< clipper::datatypes::F_phi<float> >	fphis(	original_fphis.spacegroup(),
									original_fphis.cell(),
									original_fphis.hkl_sampling()	);
	fphis = original_fphis;
	for ( hri = fphis.first(); !hri.last(); hri.next() ) {
		n_data++;
		float f = fphis[hri].f();
		if (! clipper::Util::is_nan(f)) {
			float irs =  hri.invresolsq();
			fphis[hri].f() *= exp(-b_factor * irs * 0.25); 
		}
		else {
			n_inf++;
		}
	}
	if( n_inf == n_data )
		std::cout << "INFO::ERROR" << std::endl; 
	xmap.fft_from(fphis);

	return xmap;

}

clipper::Xmap<float>
patterson( clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis ) {
	int n_inf=0,n_data=0,n_own=0;
	clipper::Grid_sampling gs( original_fphis.spacegroup(), original_fphis.cell(), original_fphis.resolution(), 1.5 );
	clipper::Xmap<float> xmap;	  					// NULL MAP
	xmap.init ( original_fphis.spacegroup(), original_fphis.cell(), gs ); 	// INITIALIZE MAP
	clipper::HKL_info::HKL_reference_index hri;
	clipper::HKL_data< clipper::datatypes::F_phi<float> >	fphis(	original_fphis.spacegroup(),
									original_fphis.cell(),
									original_fphis.hkl_sampling()	);
	fphis = original_fphis;
	for ( hri = fphis.first(); !hri.last(); hri.next() ) {
		n_data++;
		float f 	= fphis[hri].f();
		float irs 	= hri.invresolsq();
		if ( !clipper::Util::is_nan(f) ) {
			std::cout	<< "  \t  " << f << "  \t  " << irs 
					<< "\t" << hri.hkl().h() << "\t"<< hri.hkl().k() << "\t"<< hri.hkl().l() 
					<< std::endl ;n_own++;
		}
		if (! clipper::Util::is_nan(f)) {
			fphis[hri].f() = f*f; 
			fphis[hri].phi() = 0.0; 
		}
		else {
			n_inf++;
		}
	}
	if( n_inf == n_data )
		std::cout << "INFO::ERROR" << std::endl; 
	std::cout << "INFO::COUNT\t" << n_own << std::endl; 
	xmap.fft_from(fphis);
	return xmap;
}


clipper::Xmap<float>
bijvoet_diff( clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis ) {
	int n_inf=0,n_data=0,n_own=0;
	clipper::Grid_sampling gs( original_fphis.spacegroup(), original_fphis.cell(), original_fphis.resolution(), 1.5 );
	clipper::Xmap<float> xmap;	  					// NULL MAP
	xmap.init ( original_fphis.spacegroup(), original_fphis.cell(), gs ); 	// INITIALIZE MAP
	std::cout << "INFO GOT HERE! " << std::endl;
	clipper::HKL_info::HKL_reference_index hri;
	clipper::HKL_data< clipper::datatypes::F_phi<float> >	fphis(	original_fphis.spacegroup(),
									original_fphis.cell(),
									original_fphis.hkl_sampling()	);
	fphis = original_fphis;
	for ( hri = fphis.first(); !hri.last(); hri.next() ) {
		n_data++;
		float f 	= fphis[hri].f();
		float irs 	= hri.invresolsq();
		if ( !clipper::Util::is_nan(f) ) {
//
//			std::cout	<< "  \t  " << f << "  \t  " << irs 
//					<< "\t" << hri.hkl().h() << "\t"<< hri.hkl().k() << "\t"<< hri.hkl().l() 
//					<< std::endl ;n_own++;
//
			clipper::HKL rfl(hri.hkl().h()*-1,hri.hkl().k()*-1,hri.hkl().l()*-1);
			clipper::HKL mno;
			int cidx;
			//cidx	= hri.base_hkl_info().index_of(rfl);
			int sym=0; bool friedel=false; 
			mno	= hri.base_hkl_info().find_sym(rfl,sym,friedel); //index_of(rfl); 
			cidx	= hri.base_hkl_info().index_of(mno);
			if( cidx<0 ) {
				// std::cout << " " << cidx << std::endl;
				continue;
			}
			clipper::HKL_info::HKL_reference_index hri_m( hri.base_hkl_info(), cidx );
			float df	= fphis[hri_m].f();
			std::cout << "INFO \t " << df << std::endl;
			df*=df;f*=f;
			fphis[hri].f() = f-df; 
			fphis[hri].phi() = 0.0; 
		}
		else {
			n_inf++;
		}
	}
	if( n_inf == n_data )
		std::cout << "INFO::ERROR" << std::endl; 
	std::cout << "INFO::COUNT\t" << n_own << std::endl; 
	xmap.fft_from(fphis);
	std::cout << "INFO::DONE HERE" << std::endl;
 
	return xmap;
}


float optimal_B_kurtosis( clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis ) {
// 
// CALCULATES THE OPTIMAL BFACTOR:
// PERFORMS A GOLDEN RATIO SEARCH ON
// THE KURTOSIS OF THE ENTIRE SUPPLIED MAP
// 
	float sharpening_limit = 100;	
	float kurtosis = 0.0, B_optimal = 0.0;
	float a =-1.0*sharpening_limit, b=1.0*sharpening_limit, TOL=1E-2, TOLB=40E0;
	float fa= 0.0 , fb = 0.0;
	float fc= 0.0 , fd = 0.0, golden_ratio = (sqrt(5.0)-1.0)*0.5;
	float c = b-golden_ratio*(b-a);
	float d = a+golden_ratio*(b-a);
	float k = 1E-20, m=1.0;
	float a0=a;

	std::cout << "GOLDEN RATIO SEARCH::";

	if ( true ) {
		fa	= calc_kurtosis( sharpen ( a, original_fphis ), 0 );
        	fb	= calc_kurtosis( sharpen ( b, original_fphis ), 0 );
        	k = (fb-fa)/(b-a);
        	m = fa;
		while( d-c>TOL ) {
			fc		= calc_kurtosis( sharpen ( c, original_fphis ), 0 );
			fd		= calc_kurtosis( sharpen ( d, original_fphis ), 0 );
			if( d-c>TOLB ) {
				fc	/= (k*(c-a0)+m);
				fd	/= (k*(d-a0)+m);
			}
			if( fc > fd ) { // FIND MAXIMUM
				b = d; d = c;
				c = b - golden_ratio*( b - a );
			} else {
				a = c; c = d;
				d = a + golden_ratio*( b - a );
			}
		}
		B_optimal = (c+d)*0.5;
	}
	return B_optimal;
}

float optimal_B_skewness( clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis ) {
// 
// CALCULATES THE OPTIMAL BFACTOR:
// PERFORMS A GOLDEN RATIO SEARCH ON
// THE KURTOSIS OF THE ENTIRE SUPPLIED MAP
// 
	float sharpening_limit = 100;	
	float kurtosis=0.0, B_optimal=0.0;
	float a =-1.0*sharpening_limit, b=1.0*sharpening_limit, TOL=1E-2, TOLB=40E0;
	float fa= 0.0 , fb = 0.0;
	float fc= 0.0 , fd = 0.0, golden_ratio = (sqrt(5.0)-1.0)*0.5;
	float c = b-golden_ratio*(b-a);
	float d = a+golden_ratio*(b-a);
	float k = 1E-20, m=1.0;
	float a0=a;

	std::cout << "GOLDEN RATIO SEARCH::";

	if ( true ) {
		fa	= calc_skewness( sharpen ( a, original_fphis ), 0 );
        	fb	= calc_skewness( sharpen ( b, original_fphis ), 0 );
        	k	= (fb-fa)/(b-a);
        	m	=  fa;
		while( d-c > TOL ) {
			fc		= calc_skewness( sharpen ( c, original_fphis ), 0 );
			fd		= calc_skewness( sharpen ( d, original_fphis ), 0 );
			if( d-c>TOLB ) {
				fc	/= (k*(c-a0)+m);
				fd	/= (k*(d-a0)+m);
			}
			if( fc > fd ) { // FIND MAXIMUM
				b = d; d = c;
				c = b - golden_ratio*( b - a );
			} else {
				a = c; c = d;
				d = a + golden_ratio*( b - a );
			}
		}
		B_optimal = (c+d)*0.5;
	}
	return B_optimal;
}


float optimal_B_std_moment( clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis , int imom ) {
// 
// CALCULATES THE OPTIMAL BFACTOR:
// PERFORMS A GOLDEN RATIO SEARCH ON
// THE KURTOSIS OF THE ENTIRE SUPPLIED MAP
// 
	float sharpening_limit = 100;	
	float kurtosis=0.0, B_optimal=0.0;
	float a =-1.0*sharpening_limit, b=1.0*sharpening_limit, TOL=1E-2, TOLB=40E0;
	float fa= 0.0 , fb = 0.0;
	float fc= 0.0 , fd = 0.0, golden_ratio = (sqrt(5.0)-1.0)*0.5;
	float c = b-golden_ratio*(b-a);
	float d = a+golden_ratio*(b-a);
	float k = 1E-20, m=1.0;
	float a0=a;

	std::cout << "GOLDEN RATIO SEARCH::";

	if ( true ) {
		fa	= calc_standardized_moment( sharpen ( a, original_fphis ), imom, 0 );
        	fb	= calc_standardized_moment( sharpen ( b, original_fphis ), imom, 0 );
        	k = (fb-fa)/(b-a);
        	m = fa;
		while( d-c > TOL ) {
			fc		= calc_standardized_moment( sharpen ( c, original_fphis ), imom , 0 );
			fd		= calc_standardized_moment( sharpen ( d, original_fphis ), imom , 0 );
			if(d-c>TOLB) {
				fc	/= (k*(c-a0)+m);
				fd	/= (k*(d-a0)+m);
			}
			if( fc > fd ) { // FIND MAXIMUM
				b = d; d = c;
				c = b - golden_ratio*( b - a );
			} else {
				a = c; c = d;
				d = a + golden_ratio*( b - a );
			}
		}
		B_optimal = (c+d)*0.5;
	}
	return B_optimal;
}

void scale_xmap( clipper::Xmap<float>& xmap , float scalef ) {
	for( clipper::Xmap<float>::Map_reference_index midx = xmap.first(); !midx.last(); midx.next() ) {
		float x		= xmap[midx];					// THE VALUE
		xmap[midx]	= x * scalef;
	}
}

void
read_mtz_file(	std::string filename , clipper::HKL_info *myhkl ,
		clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata	,
		clipper::MTZcrystal *xtal , std::vector< std::string > labels	,
		bool match_types_not_labels )
{
	clipper::CCP4MTZfile mtzin;
	std::string	ext("mtz")	;
	std::string 	fphi_str	;
	int	is_mtz_file	= 0 ,
		verbose		= 1 ,
		write_mtz	= 0 ;

	if(verbose) {
		std::cout << "INFO::WORKING WITH " << filename <<" AND LABELS::"<< std::endl; 
		for(int i=0;i<labels.size();i++)
			std::cout << labels[i] << " " ;
		std::cout << std::endl;
	}
	std::size_t found	= filename.find(ext);

	if ( found!=std::string::npos ) 
	{
		std::cout << "INFO:: FOUND MTZ INPUT NAME" << std::endl;
		try { 
			mtzin.open_read( filename );	
			is_mtz_file = 1;
			std::cout << "INFO::OPENED FILE"<<std::endl;
		}
   		catch ( ... ) {
     			std::cout << "INFO:: NOT A VALID MTZ FILE: " << filename  << std::endl;
      			is_mtz_file = 0;
   		}

    		std::string label	;
   		std::string type	;
		std::string mtypF("F")	, mtypP("P");
		std::string mlab("WT")	;
		std::string selection( "/[" );
		std::string base_str;
		std::vector< int > found_type;
		for(int i=0;i<labels.size();i++)
			found_type.push_back(0);

		if ( is_mtz_file ) {
			int nHit = 0;
			mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
			std::vector<clipper::String> v = mtzin.column_labels();
			if ( v.size() > 1 ) {
				int read_success = 1;
				for ( unsigned int i=0 ; i<v.size() ; i++ ) {
					if(verbose)
    						std::cout << i << " " << v[i] << "\n";
   					std::string::size_type ispace  = v[i].find_last_of(" ");
   					if (ispace == std::string::npos) {
						std::cout	<< "WARNING:: uninterprettable label \" "
								<< v[i] << "\" of "<< filename[0] << "\n";
					} else {
						label	= v[i].substr ( 0 , ispace  );
						type	= v[i].substr ( ispace + 1  );
						for(int t=0 ; t<labels.size() ; t++) {
							std::size_t found ;
							if( match_types_not_labels )
								 found = type.find(	labels[t]	);
							else
								 found = label.find(	labels[t]	);

							if ( found != std::string::npos && 
								found_type[ t ] == 0	) {
								nHit++;
								base_str	 = SplitString( label , "\\/" , 0 );
								std::string	toplabel_str = SplitString( label, "\\/" , 1);
								selection	+= toplabel_str;
								selection	+= (nHit==1)?(","):("]");
								found_type[ t ]  = 1;
							}
						}
						if( verbose == 2 ) {
							std::cout << "Got label :" << label
							<< ": and type :" << type << ":\n";
							SplitFilename (label);
						}
					}
				}
			}
		}
//	HERE WE RETRIEVE DATA
		fphi_str = base_str + selection;
		// fphi_str = cols;

		if(verbose)
			std::cout << "INFO:: WILL GET DATA FROM: " << fphi_str << std::endl;
		mtzin.import_hkl_info(	( *myhkl ) );
		mtzin.import_crystal (	( *xtal  )	, fphi_str );
		mtzin.import_hkl_data(	( *fphidata )	, fphi_str );
		mtzin.close_read	( );

		if(verbose) {
			int n_reflections = fphidata->num_obs();
	      		std::cout << "INFO:: NUMBER OF OBSERVED REFLECTIONS: " << n_reflections << "\n";
			if (n_reflections <= 0 ) {
				std::cout << "WARNING:: NO REFLECTIONS IN MTZ!?" << std::endl;
				die("NOTHING TO DO");
			}
		}
	}
	std::cout << "INFO::FOUND GOOD DATA" << std::endl;
}

int main (int argc, char **argv) 
{
	int verbose = 0;
	std::string 		filename[3];
	clipper::HKL_info 	myhkl;
	clipper::HKL_data< clipper::datatypes::F_phi<float> > fphidata;
	clipper::CCP4MTZfile 	mtzin;
	clipper::MTZcrystal  	xtal;

	double sigma3d	=  5.0;
	double sigma1d	=  1.0;
	double sigma4d	= -1.0;

	std::string cols;
	std::string::size_type sz;

	std::pair < int, char ** >	margs;
	margs.first = argc; margs.second = argv;
	std::vector < std::pair < int, std::pair<std::string, std::string > > >	opts_n_defaults;
	std::pair	< int, std::pair<std::string,std::string > > 		vipss;
	vipss.first=1;	vipss.second.first	= "-i";
			vipss.second.second	= "i_map.mtz";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-o";
			vipss.second.second	= "o_map.mtz";
        opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-s1";
			vipss.second.second	= "1.0";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-s2";
			vipss.second.second	= "5.0";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-s3";
			vipss.second.second	= "10.0";
	opts_n_defaults.push_back(vipss);

	if(argc<2)
		fatal("NOT ENOUGH ARGUMENTS");
	std::vector<int> v_set;
	for( int i=0 ; i < opts_n_defaults.size() ; i++ )
		v_set.push_back(0);
	std::pair< int , std::vector< std::string > > run_args = argparser( margs, v_set, opts_n_defaults );
	if( run_args.first ) {
		fatal();
	}

	if( v_set[2] )
		sigma1d = atof( run_args.second[2].c_str() );
	if( v_set[3] )
		sigma3d = atof( run_args.second[3].c_str() );
	if( v_set[4] )
		sigma4d = atof( run_args.second[4].c_str() );

	std::string tmpline( run_args.second[0] );
	std::vector< std::string > label_types ;
	label_types.push_back("F");
	label_types.push_back("P");

	read_mtz_file( tmpline , &myhkl, &fphidata , &xtal , label_types , true );

////	DEBUG
	if( verbose == 2 ) {
		std::cout << "INFO:: MTZ:: DEBUG::    (M)YHKL:: " << myhkl.spacegroup().symbol_hm() << " "
		   << myhkl.cell().descr().a() << " "
		   << myhkl.cell().descr().b() << " "
		   << myhkl.cell().descr().c() << " "
		   << clipper::Util::rad2d(myhkl.cell().descr().alpha()) << " "
		   << clipper::Util::rad2d(myhkl.cell().descr().beta ()) << " "
		   << clipper::Util::rad2d(myhkl.cell().descr().gamma()) << " "
		   << std::endl;
		std::cout << "INFO:: MTZ:: DEBUG:: (FP)HIDATA:: " << fphidata.spacegroup().symbol_hm() << " "
		   << fphidata.cell().descr().a() << " "
		   << fphidata.cell().descr().b() << " "
		   << fphidata.cell().descr().c() << " "
		   << clipper::Util::rad2d(fphidata.cell().descr().alpha()) << " " 
		   << clipper::Util::rad2d(fphidata.cell().descr().beta ()) << " "
		   << clipper::Util::rad2d(fphidata.cell().descr().gamma()) << " "
		   << std::endl;
	}

////	THE ACTUAL PROGRAM
	const float map_sampling_rate = 1.5;
	clipper::Grid_sampling gshkl( myhkl.spacegroup(),    myhkl.cell(),    myhkl.resolution(), map_sampling_rate );  // define grid
	clipper::Grid_sampling gs( fphidata.spacegroup(), fphidata.cell(), fphidata.resolution(), map_sampling_rate );
	if ( verbose ) {
		std::cout << "INFO:: GRID SAMPLING (FP): " << gs.format()	<< std::endl; 
		std::cout << "INFO:: GRID SAMPLING (M):: " << gshkl.format()	<< std::endl; 
	}
	clipper::Xmap<float> xmap, xmap_q, xmap_p, xmap_s, xmap_m;	  	// NULL MAP
	xmap.init	( fphidata.spacegroup(), fphidata.cell(), gs); 	// INITIALIZE MAP
	xmap_q.init	( fphidata.spacegroup(), fphidata.cell(), gs); 	// INITIALIZE MAP
	xmap_p.init	( fphidata.spacegroup(), fphidata.cell(), gs); 	// INITIALIZE MAP
	xmap_m.init	( fphidata.spacegroup(), fphidata.cell(), gs); 	// INITIALIZE MAP
	xmap_s.init	( fphidata.spacegroup(), fphidata.cell(), gs); 	// INITIALIZE MAP
//
//	int smom	= 1;
//	int nmom	= 3;
//	std::vector< clipper::Xmap<float> > xmom;
//	for (int imom=0;imom<nmom;imom++){
//		clipper::Xmap<float> xmo;
//		xmo.init(fphidata.spacegroup(), fphidata.cell(), gs);
//		xmom.push_back(xmo);
//	}
//
////	DO FFT
	if(verbose)
		std::cout << "FFT :: BEGUN" << std::endl;
	xmap.fft_from( fphidata );
	if(verbose)                  			 	// GET ERROR FROM CLIPPER
		std::cout << "FFT :: ENDED" << std::endl;
	float fnorm	= 0.0;
	float fnorm0	= 0.0;
	std::pair<float,float> mean_std0, mean_std1, mean_std2;

	// for( int i=0 ; i<-1 ; i++ )
	{
		int i	 = 0 ;
		float dc = 0.40 + ((float) i ) * 0.005;
		for( clipper::Xmap<float>::Map_reference_index midx = xmap.first(); !midx.last(); midx.next() ){
			xmap_q[midx]=xmap[midx]<dc?xmap[midx]:0.0;
		}
		std::cout << " @ \t  " << dc << " \t "
			<< calc_standardized_moment( xmap_q , 1 , 0 ) << " \t "
			<< calc_standardized_moment( xmap_q , 2 , 0 ) << " \t "
                        << calc_standardized_moment( xmap_q , 3 , 0 ) << " \t "
			<< "\n" ;
	}
//
//	xmap_m += xmap;
//	mean_std0 = map_mean_variance(xmap,0);
//	fnorm0 = sum_map( xmap_m , 0 );
//
	std::vector< float > stg;
	float bopt = optimal_B_kurtosis(	fphidata 	);
	xmap_s	= sharpen	( 	bopt,	fphidata 	);
	xmap_q	= degaus	( 	xmap_s , stg 		);
//
//	} else {
//		float bopt_k	= optimal_B_kurtosis( fphidata );
//		xmap_q		= sharpen ( bopt_k ,  fphidata );
//		float bopt_s	= optimal_B_skewness( fphidata );
//		xmap_s		= sharpen ( bopt_s ,  fphidata );
//		std::cout << "\n INFO [ " << bopt_s << " \t " << bopt_k << " ] " << std::endl;
//		float bopt[nmom];
//		std::cout	<< " SMOM [ " ;
//		for( int imom=smom; imom<nmom; imom++ ) {
//			bopt[imom]	= optimal_B_std_moment( fphidata, imom );
//			xmom[imom]	= sharpen ( bopt[imom], fphidata );
//			if( 1 ) {
//				xmap_m	+=  xmom[imom];
//			}
//			std::cout << " \t " << imom << " " << bopt[imom] << "\n";
//		}
//		std::cout << " ]" <<std::endl;
//	}
//	fnorm = sum_map( xmap_m , 0 ); 
//	float norm=(float)nmom;
//	mean_std1=map_mean_variance(xmap_m,0);

//// RENORMALIZATION
//	scale_xmap( xmap_m , mean_std0.second/mean_std1.second );
//	mean_std2 	= map_mean_variance(xmap_m,0);
//	xmap_p 		= degaus( xmap_m, 5.0 , mean_std2.second , false );
	std::string base_str = SplitString( filename[0] , "." , 0 );

////	FILE IO
	clipper::CCP4MAPfile mapout;
	filename[1]  =  tmpline  + "_bk.map";
	mapout.open_write	( filename[1] 	);
	mapout.export_xmap	( xmap_q );
	mapout.close_write	(	 );

	filename[1]  =  tmpline  + "_0.map";
	mapout.open_write	( filename[1] 	);
	mapout.export_xmap	( xmap );
	mapout.close_write	(	);

	if( !( sigma1d>0.0 && sigma3d>0.0 ) ) {
		filename[1] = base_str + "_bs.map";
		mapout.open_write	( filename[1] 	);
		mapout.export_xmap	( xmap_p );
		mapout.close_write	(	 );
		clipper::CCP4MAPfile mapout;
		filename[1] = base_str + "_skew.map";
		mapout.open_write	( filename[1] 	);
		mapout.export_xmap	( xmap_q );
		mapout.close_write	(	 );
		filename[1] = base_str + "_ms.map";
		mapout.open_write	( filename[1] 	);
		mapout.export_xmap	( xmap_m );
		mapout.close_write	(	 );
	}

	return 0;
}
