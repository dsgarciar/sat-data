
struct {	char	line1	[70];
			char	line2	[70];
			char	name	[25];
			long	catnum;
			long	setnum;
			char	designator	[10];
			int		year;
			double	incl;
			double	raan;
			double	eccn;
			double	argper;
			double	meanan;
			double	meanmo;
			double	drag;
			double	nddot6;
			double	bstar;
			long	orbitnum;
	}	sat	[24];

struct	{	char	callsign	[17];
			double	stnlat;
			double	stnlong;
			int		stnalt;
	}	qth;

struct	{	char	name	[25];
			long	catnum;
			char	squintflag;
			double	alat;
			double	alon;
			unsigned	char	transponders;
			char	transponder_name	[10][80];
			double	uplink_start	[10];
			double	uplink_end	[10];
			double	downlink_start	[10];
			double	downlink_end	[10];
			unsigned	char	dayofweek	[10];
			int		phase_start	[10];
			int		phase_end	[10];
	}	sat_db	[24];

double	tsince, jul_epoc, jul_utc, eclipse_depth=0,
		sat_azi, sat_ele, sat_range, sat_range_rate,
		sat_lat, sat_lon, sat_alt, sat_vel, phase,
		sun_azi, sun_ele, daynum, fm, fk, age, aostime,
		lostime, ax, ay, az, rx, ry, rz, squint, alat, alon,
		sun_ra, sun_dec, sun_lat, sun_lon, sun_range, sun_range_rate,
		moon_az, moon_el, moon_dx, moon_ra, moon_dec, moon_qha, moon_dv;

char	qthfile[50], tlefile[50], dbfile[50], temp[80], output[25],
		serial_port[15], resave=0, reload_tle=0, netport[7],
		once_per_second=0, ephem[5], sat_sun_status, findsun,
		calc_squint, database=0, xterm, io_lat='N', io_lon='W';

int		indx, antfd, iaz, iel, ma256, isplat, isplong, socket_flag=0,
		Flags=0;

long	rv, irk;

unsigned	char	val	[256];

char	visibility_array[24], tracking_mode[30];

float	az_array[24], el_array[24], long_array[24], lat_array[24],
		footprint_array[24], range_array[24], altitude_array[24],
		velocity_array[24], eclipse_depth_array[24], phase_array[24],
		squint_array[24];

double	doppler[24], orbitnum_array[24];

unsigned	short	portbase=0;

typedef	struct	{
		double	epoch, xndt2o, xndd6o, bstar, xincl,
			xnodeo, eo, omegao, xmo, xno;
		int		catnr, elset, revnum;
		char	sat_name[25], idesg[9];
	} tle_t;

typedef	struct	{
		doube lat, lon, alt, theta;
	}	geodetic_t;

typedef	struct	{
		double x, y, z, w;
	}	vector_t;

typedef	struct	{
		double	eosq, sinio, cosio, betao, aodp, theta2,
			sing, cosg, betao2, xmdot, omgdot, xnodot, xnodp;
		double	x11, omgadf, xnode, em, xinc, xn, t;
		double	ds50;
	}	deep_arg_t;

geodetic_t	obs_geodetic;

tle_t	tle;


