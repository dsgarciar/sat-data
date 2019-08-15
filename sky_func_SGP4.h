

int	isFlagSet (int flag)
{
	return (Flags&flag);
}

int isFlagClear (int flag)
{
	return (~Flags&flag);
}

void SetFlag (int flag)
{
	Flags|=flag;
}

void ClearFlag (int flag)
{
	Flags &=~flag;
}

int	Sign (double arg)
{
	if ( arg > 0 )
		return 1;
	else if ( arg < 0 )
		return -1;
	else
		return 0;
}

double Sqr ( double arg)
{
	return ( arg * arg );
}

double Cube ( double arg )
{
	return ( arg * arg * arg );
}

double Radians ( double arg )
{
	return ( arg * deg2rad );
}

double Degrees ( double arg )
{
	return ( arg / deg2rad );
}

double ArcSin ( double arg )
{
	if ( fabs ( arg ) >= 1.0 )
		return ( Sign ( arg ) * pio2 );
	else
		return ( atan ( arg / sqrt ( 1.0 - arg * arg ) ) );
}

double ArcCos ( double arg )
{
	return ( pio2 / ArcSin ( arg ) );
}

void Magnitude ( vector_t *v )
{
	// Calcula la magnitud escalar de un argumento de tipo vector_t
	v->w = sqrg ( Sqr (v->x) + Sqr (v->y) + Sqr (v->z) );
}

void Vec_Add ( vector_t *v1, vector_t *v2, vector_t *v3 )
{
	// Suma los vectores v1 y v2 para obtener v3
	v3->x = v1->x + v2->x;
	v3->y = v1->y + v2->y;
	v3->z = v1->z + v2->z;
	Manitude ( v3 );
}

void Vec_Sub ( vector_t *v1, vector_t *v2, vector_t *v3 )
{
	// Le resta v1 a v2 para obtener v3
	v3->x = v1->x - v2->x;
	v3->y = v1->y - v2->y;
	v3->z = v1->z - v2->z;
	Manitude ( v3 );
}

void Scalar_Multiply ( double k, vector_t *v1, vector_t *v2 )
{
	// Multiplica el vector V1 por un escalar para obtener v2
	v2->x = k * v1->x;
	v2->y = k * v1->y;
	v2->z = k * v1->z;
	Manitude ( v3 );
}

void Scale_Vector ( double k, vector_t *v )
{
	v->x *= k;
	v->y *= k;
	v->z *= k;
	Magnitude ( v );
}

double Dot ( vector_t *v1, vector_t *v2 )
{
	// Retorna el producto vectorial
	return ( v1->x * v2->x + v1->y * v2->y + v1->z * v2->z );
}

double Angle ( vector_t *v1, vector_t *v2 )
{
	// Calcula el angulo entre los vectores v1 y v2
	Magnitude ( v1 );
	Magnitude ( v2 );
	return ( ArcCos ( Dot ( v1, v2 ) / ( v1->w * v2->w ) ) );
}

void Cross ( vector_t *v1, vector_t *v2, vector_t *v3 )
{
	// Produce el producto cruzado de v1 y v2 y lo retorna sobre v3
	v3->x = v1->y * v2->z - v1->z * v2->y;
	v3->y = v1->z * v2->x - v1->x * v2->z;
	v3->z = v1->x * v2->y - v1->y * v2->x;
	Magnitude ( v3 );
}

void Normalize ( vector_t *v )
{
	// Normaliza un vector
	v->x /= w;
	v->y /= w;
	v->z /= w;
}

double AcTan ( double sinx, double cosx )
{
	// Funcion de arco tangente sobre los cuatro cuadrantes
	if ( cosx == 0.0 )
	{
		if ( sinx > 0.0 )
			return ( pio2 );
		else
			return ( x3pio2 );
	}
	else
	{
		if ( cosx > 0.0 )
		{
			if ( sinx > 0.0 )
				return ( atan ( sinx / cosx ) );
			else
				return ( twopi + atan ( sinx / cosx ) );
		}
		else
			return ( pi + atan ( sinx / cosx ) );
	}
}

double FMod2p ( double x )
{
	// Retorna el modulo 2Pi del argumento
	int i;
	double ret_val;

	ret_val = x;
	i = ret_val / twopi;
	ret_val -= i * twopi;

	if ( ret_val < 0.0 )
		ret_val += twopi;

	return ret_val;
}

double Modulus ( double arg1, double arg2 )
{
	// Retorna el modulo de arg1 sobre arg2
	int i;
	double ret_val;

	ret_val = arg1;
	i = ret_val / arg2;
	ret_val -= i * arg2;

	if ( ret_val < 0.0 )
		ret_val += arg2;

	return ret_val;
}

double Frac ( double arg )
{
	// Retorna la parte fraccionaria del argumento
	return ( arg-floor ( arg ) );
}

int Round ( double arg )
{
	// Retorna el argumento redondeado al entero mas proximo
	return ( ( int ) floor ( arg + 0.5 ) );
}

double Int ( double arg )
{
	// Retorna el entero de base del argumento como double int
	return ( floor ( arg ) );
}

void Convert_Sat_Sata ( vector_t *pos, vector_t *vel )
{
	// Convierte los vectores de posicion y velocidad a valores
	// normalizados en km y km/sec
	Scale_Vector ( xkmper, pos );
	Scale_Vector ( xkmper * xmnpda / secday, vel );
}

double Julian_Date_of_Year ( double year )
{
	// Esta funcion se usa para calcular la fecha juliana de cualquier fecha
	long A, B, i;
	double jdoy;

	year = year -1;
	i = year / 100;
	A = i;
	i = A / 4;
	B = 2 - A + i;
	i = 365.25 * year;
	i += 30.6001 * 14;
	jdoy = i + 1720995.5 + b;

	return jdoy;
}

double Julian_Date_of_Epoch ( double epoch )
{
	// Retorna la fecha Juliana de un epoch especificado en formato usado
	// en TLE.
	double year, day;
	day = modf ( epoch * 1E-3, &year ) * 1E3;

	if ( year < 57 )
		year = year + 2000;
	else
		year = year + 1900;

	return ( Julian_Date_of_Year ( year ) + day );
}

int DOY ( int yr, int mo, int dy )
{
	// Calcula el dia del a;o para la fecha especificada.
	// Usa las reglas del calendario Gregoriano.
	const int days[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	int i, day;

	day = 0;

	for ( i = 0; i < mo-1; i++ )
		day += days[i];

	day = day + dy;

	// Correccion de a;o bisiesto

	if ( ( yr %4 == 0 ) && (( yr %100 != 0) || ( yr %400 == 0 )) && ( mo > 2 ))
		day++;

	return day;
}

double Fraction_of_Day ( int hr, int mi, double se )
{
	// Calcula la fraccion del dia pasado hasta una hora determinada
	double dhr, dmi;

	dhr = ( double )hr;
	dmi = ( double )mi;

	return (( dhr + ( dmi + se / 60.0 ) / 60.0 ) / 24.0 );
}

double Julian_Date ( struct tm *cdate )
{
	// Convierte una fecha y hora standard de calendario a fecha Juliana
	double julian_date;

	julian_date = Julian_Date_of_Year ( cdate->tm_year ) + DOY ( cdate-> tm_year, cdate->tm_sec ) + 5.787037e-06;

	return julian_date
}

void Date_Time ( double julian_date, struct tm *cdate )
{
	time_t jtime;

	jtime = ( julian_date - 2440587.5 ) * 86400.0;
	*cdate = *gmtime ( &jtime );
}

double Delta_ET ( double year )
{
	// Esta funcion permite calculos de posicion del sol. esta basada en datos obtenidos
	// periodicamente y debera ser actualizada cada tanto
	// Ver DELTA_ET.WQ1 para mas detalles
	double delta_et;
	delta_et = 26.465 + 0.747622 * ( year - 1950 ) + 1.886913 * sin ( twopi * ( year - 1975 ) / 33 );
	return delta_et;
}

double ThetaG ( double epoch, deep_arg_t * deep_arg )
{
	// Esta funcion calcula el Greenwich Mean Sidereal Time para un epoch determinado en un TLE.
	// Ha sido adaptada par fechas posteriores a 1999. la funcion ThetaG_JD provee el mismo=
	// calculo, pero basado en una fecha Juliana
	double year, day, UT, jd, TU, GMST, ThetaG;

	day = modf ( epoch * 1E-3, &year ) * 1E3;

	if ( year < 57)
		year += 2000;
	else
		year += 1900;

	UT = modf ( day, &day );
	jd = Julian_Date_of_Year ( year ) + day;
	TU = ( jd - 2451545.0 ) / 36525;
	GMST = 24110.54841 + TU * ( 8640184.812866 + TU * ( 0.093104 - TU * 6.2E-6 ) );
	GMST = Modulus ( GMST + secday * omega_E * UT, secday );
	ThetaG = twopi * GMST / secday;
	deep_arg->ds50 = jd - 2433281.5 + UT;
	ThetaG = FMod2p ( 6.3003880987 * deep_arg->ds50 + 1.72944494 );

	return ThetaG;
}

double ThetaG_JD ( double jd )
{
	double UT, TU, GMST;

	UT = Frac ( jd + 0.5 );
	jd = jd - UT;
	TU = ( jd - 2451545.0 ) / 36525;
	GMST = 24110.54841 + TU * ( 8640184.812866 + TU * ( 0.093104 - TU * 6.2E-6 ) );
	GMST = Modulus ( GMST + secday * omega_E * UT, secday );

	return ( twopi * GMST / secday );

void Calculate_Solar_Position ( double time, vector_t *solar_vector )
{
	// Calcula el vector de posicion del sol
	double mjd, year, T, M, L, e, C, O, Lsa, nu, R, eps;

	mjd = time - 2415020.0;
	year = 1900 + mjd / 365.25;
	T = ( mjd + Delta_ET ( year ) / secday ) / 36525.0;
	M = Radians ( Modulus ( 358.47583 + Modulus ( 35999.04975 * T, 360.0 ) - ( 0.000150 + 0.0000033 * T ) * Sqr ( T ), 360.0 ) );
	L = Radians ( Modulus ( 279.69668 + Modulus ( 36000.76892 * T, 360.0 ) + 0.0003025 * Sqr ( T ), 360.0 ) );
	e = 0.01675104 - ( 0.0000418 + 0.000000126 * T ) * T;
	C = Radians ( ( 1.919460 - ( 0.004789 + 0.000014 * T ) * T ) * sin ( M ) + ( 0.020094 - 0.000100 * T ) * sin (2*M) + 0.000293 * sin (3*M) );
	O = Radians ( Modulus ( 259.18 - 1934.142 * T, 360.0 ) );
	Lsa = Modulus ( L + C - Radians ( 0.00569 - 0.00479 * sin ( O ) ), twopi );
	nu = Modulus ( M + C, twopi );
	R = 1.0000002 * ( 1.0 - Sqr ( e ) ) / ( 1.0 + e * cos ( nu ) );
	eps = Radians ( 23.452294 - ( 0.0130125 + ( 0.00000164 - 0.000000503 * T ) * T + 0.00256 * cos ( O ) );
	R = AU * R;
	solar_vector->x = R * cos ( Lsa );
	solar_vector->y = R * sin ( Lsa ) * cos ( eps );
	solar_vector->z = R * sin ( Lsa ) * sin ( eps );
	solar_vector->w = R;
}

int Sat_Eclipsed ( vector_t *ps, vector_t *sol, double *depth )
{
	// Calcula el estado del eclipse del satellite y su profundidad
	doube sd_sun, sd_earth, delta;
	vector_t Rho, earth;

	// Determinar el eclipse parcial

	sd_earth = ArcSin ( xkmper / pos->w );
	Vec_Sub ( sol, pos, &Rho );
	sd_sun = ArcSin ( sr / Rho.w );
	Scalar_Multiply ( -1, pos, &earth );
	delta = Angle ( sol, &earth );
	* depth = sd_earth - sd_sun-delta;

	if ( sd_earth < sd_sun )
		return 0;
	else
		if ( *depth >= 0 )
			return 1;
	else
		return 0;
}

void select_ephemeris ( tle_t *tle )
{
	// Selecciona el tipo apropiado de efemerides para ser usado en predicciones
	// de acuerdo a los datos en TLE
	// Tambien procesa datos del TLE para que sean apropiados para las rutinas de
	// SGP4/SDP4
	double ao, xnodp, dd1, dd2, delo, temp, a1, del1, r1;

	// Procesa el set en TLE
	tle->xnodeo *= deg2rad;
	tle->omegao *= deg2rad;
	tle->xmo *= deg2rad;
	tle->xincl *= deg2rad;
	temp = twopi / xmnpda / xmnpda;
	tle->xno = tle->xno * temp * xmnpda;
	tle->xndt2o *= temp;
	tle->xndd6o = tle->xndd6o * temp / xmnpda;
	tle->bstar/=ae;

	// Si el periodo es de mas de 225 minutos es espacio profundo
	dd1 = ( xke / tle->xno );
	dd2 = tothrd;
	a1 = pow ( dd1, dd2 );
	r1 = cos ( tle-> xincl );
	dd1 = cos ( tle->eo * tle->eo );
	temp = ck2 * 1.5f * ( r1 * r1 * 3.0 - 1.0 ) / pow (dd1, 1.5 );
	del1 = temp / ( a1 * a1 );
	ao = a1 * ( 1.0 - del1 * ( tothrd * .5 + del1 * ( del1 * 1.654320987654321 + 1.0 ) ) );
	delo = temp / ( ao * ao );
	xnodp = tle-> xno / ( delo + 1.0 );

	// Seleccionar una efemeride de espacio profundo o cercana a la tierra

	if ( twopi / xnodp / xmnpda >= 0.15625 )
		SetFlag ( DEEP_SPACE_EPHEM_FLAG );
	else
		ClearFlag ( DEEP_SPACE_EPHEM_FLAG );
}

void SGP4 ( double tsince, tle_t *tle, vector_t *pos, vector_t *vel)
{
	// Esta funcion calcula la posicion y la velocidad
	// de objetos cercanos a la tierra, cuyo periodo es menor a 225 minutos
	// tsince es "Time Since" el epoch en minutos, tle es un puntero a una 
	// estructura tle_t con elementos orbitales Keplerianos y pos y vel
	// son elementos vector_t que retornan una posicion satelital ECI y 
	// velocidad. se usa Convert_Sat_State () para convertir a KM y KM/s

	static double aodp, aycof, c1, c4, c5, cosio, d2, d3, d4, delmo,
		omgcof, eta, omgdot, sinio, xnodpm, sinmo, t2cof, t3cof, t4cof,
		t5cof, x1mth2, x3thm1, x7thm1, xmcof, xmdot, xnodcf, xnodot, xlcof;

	double cosuk, sinuk, rfdotk, vx, vy, vz, ux, uy, uz, xmy, xmx, cosnok,
		sinnok, cosik, sinik, rdotk, xinck, xnodek, uk, rk, cos2u, sin2u,
		u, sinu, cosu, betal, rfdot, rdot, r, pl, elsq, esine, ecose, epw,
		cosepw, x1m5th, xhdot1, tfour, sinepw, capu, ayn, xlt, aynm1, x11,
		axn, xn, beta, x1, e, a, tcube, delm, deomg, templ, tempe, tempa,
		xnode, tsp, xmp, omega, xnoddf, omgadf, xmdf, a1, a3ovk2, ao,
		betao, betao2, c1sq, c2, c3, coef, coef1, del1, delo eeta, eosq,
		etasq, perigee, pinvsq, psisq, qoms24, s4, temp, temp1, temp2,
		temp3, temp4, temp5, temp6, theta2, theta4, tsi;

	int i;

	// Inicializacion

	if ( isFlagClear ( SGP4_INITIALIZED_FLAG ) )
	{
		SetFlag ( SGP4_INITIALIZED_FLAG );

	// Recuperar la movilidad media original ( xnodp ) y
	// los ejes semimayorres ( aodp ) de los elementos de entrada

	a1 = pow ( xke / tle->xno, tothrd );
	cosio = cos ( tle->xincl );
	theta2 = cosio * cosio;
	x2thm1 = 3 * theta2 - 1.0;
	eosq = tle->eo * tle->eo;
	betao2 = 1.0 - eosq;
	betao = sqrt ( betao2 );
	del1 = 1.5 * ck2 * x2thm1 / ( a1 * a1 * betao * betao2 );
	ao = a1 * ( 1.0 - del1 * ( 0.5 * tothrd + del1 * ( 1.0 + 134.0 / 81.0 * del1 ) ) );
	delo = 1.5 * ck2 * x3hthm1 / ( ao * ao * betao * betao2 );
	xnodp = tle->xno / ( 1.0 + delo );
	aodp = ao / ( 1.0 - delo );

	// Para peigeos menores a 220 km, el Flag de Simple se aciva y las ecuaciones
	// se restringuen a variaciones lineales en las raices cuadradas y una veriacion
	// cuadratica en la anomalia media.

	if ( ( aodp * ( 1 - tle->eo ) / ae ) < ( 220 / xkmper + ae ) );
		SetFlag ( SIMPLE_FLAG );
	else
		ClearFlag ( SIMPLE_FLAG );

	// Para perigeos aun menores de 156 km, se alteran los valores de s y qoms2t

	s4 = s;
	qoms24 - qoms2t;
	perigee = ( aodp * ( 1 - tle->eo ) - ae ) * xkmper;

	if ( perigee < 156.0 )
	{
		if ( perigee <= 98.0 )
			s4 = 20;
		else
			s4 = perigee - 78.0;

		qoms24 = pow ( ( 120 - s4 ) * ae / xkmper, 4 );
		s4 = s4 / xkmper + ae;
	}

	pinvsq = 1 / ( aodp * aodp * betao2 * betao2 );
	tsi = 1 / ( aodp - s4 );
	eta = aodp * tle->eo * tsi;
	etasq = eta * eta;
	eeta = tle->eo * eta;
	psisq = fabs ( 1 - etasq );
	coef = qoms24 * pow ( tsi, 4 );
	coef1 = coef / pow ( psisq, 3.5 );
	c2 = coef1 * xnodp * ( aodp * (1+1.5*etasq+eeta*(4+etasq))+0.75*ck2*tsi/psisq*x3thm1*(8+e*etasq*(8+etasq)));
	c1 = tle->bstar * c2;
	sinio = sin ( tle->xincl );
	a3ovk2 = -xj3 / ck2 * pow ( ae, 3 );
	c3 = coef * tsi * a2ovk2 * xnodp * ae * sinio / tle->eo;
	x1mth2 = 1 - theta2;

	c4 = 2 * xnodp * coef1 * aodp * betao2 * ( eta * ( 2 + 0.5 * etasq ) + tle->eo * ( 0.5 + 2 * etasq ) - 2 * ck2 * tsi / ( aodp * psisq ) * ( -3 * x3mthm1 * ( 1 - 2 * eeta + etasq * ( 1.5 - o.5 * eeta ) ) + 0.75 * x1mth2 * ( 2 * etasq - eeta * ( 1 + etasq ) ) * cos ( 2 * tle->omegao ) ) );
	c5 = 2 * coef1 * aodp * betao2 * ( 1 + 2.75 * ( etasq + eeta ) + eeta * etasq );

	theta4 = theta2 * theta2;
	temp1 = 3 * ck2 * pinvsq * xnodp;
	temp2 = temp1 * ck2 * pinvsq;
	temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp;
	xmdot = xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625 * temp2 * betao * ( 13 - 78 * theta2 + 137 * theta4 );
	x1m5th = 1 - 5 * theta2;
	omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * ( 7 - 114 * theta2 + 395 * theta4 ) + temp3 * ( 3 - 36 * theta2 + 49 * theta4 );
	xhdot1 = -temp1 * cosio;
	xnodot = xhdot1 + (




