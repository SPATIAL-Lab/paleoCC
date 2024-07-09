/*
12/11/09
-Added switches for weathering, marine OC, and biosphere feedbacks
-Made CO2 fertilization a saturating feedback

12/8/09
-Returned to perscribed warming to reduce instability
-Added Q10 for marine organic carbon burial in "warm" function

11-10-09
-Tweaked injection function
-Modified biosphere respiration feedback to be Q10 function of temperature
-Modified warm to follow CO2 sensitivity factor

10-17-09
-Fixed "rivorg" accounting - C was being added to atmosphere while isotopes
were going into the surface ocean.
-Adjusted "calburial" to reflect actual net burial of calcite and was able to
confirm steady-state mass balance within 0.02% and balance of kerogen weathering
and organic burial to within 2%.
-Version saved as GoodModern091017 has nice tuning for pCO2 of 350 ppmv.
-Need to work on re-tuning to high CO2 by changing silicate and carbonate
weathering setpoints.

10-15-09
-Found order-of-magnitude errors in atmco2 and biosphere.
-Retuned for higher silicate weathering flux (closer to original Walker values).
-Tried tuning for lower organic burial rates (ref modern marine C values of
0.17% or 5.3 x 10^12 mol for ocean productivity of 3,000 x 10^12 mol/yr;
10% of riverine POC export or 6.3 x 10^12 mol/yr; Schlunz and Schneider, 2000).
Lower orgburial (0.5%) for marine C gives shallow lysocline (1.5km).
Combining with lower river org burial (20%) gives nice case trending towards
high pCO2 (~4x) and normal lysocline (~4 km).
-Tested 1 kyr injection of 3000 GT at -40 per mil, produces excellent
2 km lysocline shoaling and 3 per mil excursion!!!!

3-17-02
Fixed dissolution flux rate problem.  Retuned for I/O ratio of 0.25.  Get
much improved C/S ratios for sediments.

2-27-2003
Updated carbon cycle using numerical recepies 6th-order Runge-Kutta solver.
In addition to cosmetic improvements, I fixed several problems in the 2001
code.  These include: 1) sediment mass transfer is based on the appropriate
equations for conservation of volume, 2) the ODEs for sediment parameters
have been corrected, no forward-looking artifical stabilization is included,
and 3) inaccuracies in the sedimentary carbon isotope equations have been
corrected.  Other major changes include the addition of several functions
for parameter calculations and the storage of basin geometry parameters in
a global matrix.  A 3-step injection and single warming are included.

6-2001
Carbon and phosphorus cycle with lysocline and variable calcite burial.
Dynamic production as a function of phosphorus concentration
Riverine fluxes of carbon and P.  Includes output to file.
Carbon isotopes.  Shelf area.  Weathering feedback and silicate weathering.
Atmosphere-surface ocean mass exchange is mass independent.
code for 15 x 5 seafloor sediment reservoirs...integrated sediments.
*/
#include "nr.h"
#include <new.h>
//#include <condefs>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#pragma hdrstop
//---------------Global Variables and Arrays---------------------------------
//USEUNIT("rkqs.cpp");
//USEUNIT("rkck.cpp");
//USEUNIT("odeint.cpp");
//---------------------------------------------------------------------------
const int sfcells = 15;      //divide ocean basins into sfcells depth intervals
const int sfdepths = 25;      //sfdepths reservoirs per depth interval
const double sfthick = 0.5 / sfdepths;
const int nvars = 12;
const int nrow = nvars + sfcells * sfdepths * 3;    //# of equations & unknowns
const int ncol = nrow + 1;
int ic = 0;
double ti, tw;

ofstream fout;
ifstream fin;

double sfres[sfcells][4];

double sco3, dco3, initsco3, initdco3, spco2, dpco2, etimesa, sigcriv, calburial, caldiss,
alkriv, fracsa, c13diss, rivcal, rivsil, rivorg, lyso, exprod, shelfc, assim, resp, eros,
fpoc, cpoc, doc, biofrac, shelfdiss, shelfdiss13C, hydro, po4riv, initCO2, initbio;
double atmco2 = 0.0495;             //10^18 mol
double vdeep = 1.23;               //10^18 m^3
double vsurf = 0.12;                //10^18 m^3
double mixing = 0.001;             //10^18 m^3/y

double shelfratio = 2.0;
double exratio = 0.040;      //fraction of standing biomass exported/year
double prodfac = 1.0;
double ioratio = 0.25;         //0.18,0.25
double particle = 0.02;        //large particle preservaton flux
double orgburial = 0.01;        //0.010, 0.005
double oxyresp = 0.01;         //fraction of respiration leading to calcite dissolution on the sea floor
double cpratio = 0.104;            //redfield c/p ratio / 1000
double fracorg = 24.0;
double fraccal = 0.0;
double fpocburial = 0.4;     //0.4...0.45,0.2, 0.4

double cwz = 1.2e-5 / 1.75;          //1.2...1.05,1.5,1.2
double swz = 0.60e-5 / 1.75;        //0.8...0.28,0.5,0.405,0.36,0.44,0.5,0.6
double owz = 0.35e-5 / 1.75;      //0.5...0.24,0.3,0.4,0.5
double pwz = 1.353e-5 / 1.25;    //1.8658...0.3,0.25,0.22,3.9,2.0,1.9,1.6,1.35
double c13rivorg = -22.5;
double c13rivcal = 2.6;

double volcano = 0.4e-5;
double c13volcano = -0.7;

double inject = 0.0;
double c13inject = -55.0;

double swtemp = 290.0;                 // initial temp in K
double dwtemp = 281.0;

//options
int inj = 1;
int fb_bio = 1;
int fb_oc = 1;
int fb_weath = 1;

Vec_DP* xp_p;
Mat_DP* yp_p;
DP dxsav;
int kmax, kount;

//---------------Function Prototypes-----------------------------------------
void initial(Vec_IO_DP& y);
void equations(const double time, Vec_I_DP& y, Vec_IO_DP& dydx);
void update(Vec_I_DP& y, const double time);
void injection(const double time);
void warm(const double time);
void carbeq(const double temp, const double sigc, const double alk,
    double& co3, double& pco2);
void lysocline();
void production(const double spo4);
void fractionation();
void biology(const double pco2, const double bioC);
void burial(Vec_I_DP& y, Vec_IO_DP& dydx);
double porosity(double c, double s, int ind);
double Ds(double vn, double c, double ct, double s, double st, double pores, double& dc, int ind);
//---------------------------------------------------------------------------
#pragma argsused
int main(int argc, char** argv)
{
    clock_t first = clock();
    Vec_DP y(nrow);
    Vec_DP dydx(nrow);
    double t1 = 0;
    double tstep = 1.0;
    double finish = 250000.0;
    double mstep = 0.001;
    dxsav = 1000.0;
    double eps = 1.0e-4;
    int nok, nbad;
    xp_p = new Vec_DP(finish / dxsav + 1);
    yp_p = new Mat_DP(nrow + 13, finish / dxsav + 1);
    kmax = finish / dxsav + 1;
    Vec_DP& xp = *xp_p;
    Mat_DP& yp = *yp_p;

    initial(y);

    int ii = 0;		//0 for simulation, 1 for spinup
    switch (ii) {

    case 0:             //for simulations
        NR::odeint(y, t1, finish, eps, tstep, mstep, nok, nbad, equations, NR::rkqs);
        break;

    case 1:             //for spinup
        for (int blech = 0; 10 - blech; blech++) {

            NR::odeint(y, t1, finish, eps, tstep, mstep, nok, nbad, equations, NR::rkqs);

            char sf[2][11] = { "tmp00.txt", "tmp01.txt" };
            for (int i = 0; 2 - i; i++) {
                fout.open(sf[i]);
                for (int j = 0; j - sfcells * 3; j++)
                    fout << yp[j + 12 + i * sfcells * 3][kount - 1] << "\n";
                fout.close();
            }
            char sedfile[10] = "tmp00.txt";
            for (int counter2 = 0; sfdepths - counter2; counter2++)
            {
                if (counter2)
                    sedfile[4] = '1';
                fin.open(sedfile);
                for (int counter = 0; sfcells - counter; counter++)
                    fin >> y[nvars + 45 * counter2 + counter];
                for (int counter = 0; sfcells - counter; counter++)
                    fin >> y[nvars + sfcells + 45 * counter2 + counter];
                fin.close();
                for (int counter = 0; sfcells - counter; counter++)
                    y[nvars + sfcells * 2 + 45 * counter2 + counter] = y[8] - fraccal;
            }

            cout << nok << "\n";
        }
        break;
    }

    //output C cycle timeseries
    fout.open("exogenic.txt");
    fout << "Time\tPCO2\tSigCSurf\tSigCDeep\tAlkSurf\tAlkDeep"
        << "\tPO4Surf\tPO4Deep\td13CAtm\td13CSurf\t"
        << "d13CDeep\tBioC\td13CBio\tInject\tStemp\tDtemp\tExport\t"
        << "Lysocline\tCalBurial\tDissolution\tc13Diss\tOrgBurial\tPO4Burial"
        << "\tsCO3--Change\tdCO3--Change\tSAFrac\n";
    for (int i = 0; kount - i; i++) {
        fout << xp[i] << "\t";
        for (int j = 0; nvars - j; j++) fout << yp[j][i] << "\t";
        for (int j = nrow; j - (nrow + 13); j++) fout << yp[j][i] << "\t";
        fout << "\n";
    }
    fout << "number of steps\t" << nok << "\n";
    fout << "elapsed time\t" << (clock() - first) / CLOCKS_PER_SEC;
    fout.close();

    //output sedimentary C timeseries
    char nfile[15] = "sediment00.txt";
    for (int j = 0; sfdepths - j; j++) {
        nfile[8] = char(j / 10 + 48);
        nfile[9] = char(j % 10 + 48);
        fout.open(nfile);
        fout << "Time\t";
        for (int i = 0; sfcells - i; i++)
            fout << "Carb" << (i + 1) << "\t";
        for (int i = 0; sfcells - i; i++)
            fout << "Sil" << (i + 1) << "\t";
        for (int i = 0; sfcells - i; i++)
            fout << "c13" << (i + 1) << "\t";
        fout << "\n";
        nfile[9] = char(j + 48);
        for (int i = 0; kount - i; i++) {
            fout << xp[i] << "\t";
            for (int k = 0; k - sfcells * 3; k++)
                fout << yp[k + nvars + j * sfcells * 3][i] << "\t";
            fout << "\n";
        }
        fout.close();
    }

    //    cout << '\a';
    return 0;
}
//---------------------------------------------------------------------------
void initial(Vec_IO_DP& y)
{
    /*     y[0] = 2.44;                 //atm pCO2
        y[1] = 2.90;               //sigC s mol/m^3
        y[2] = 3.06;                //sig C deep mol/m^3
        y[3] = 3.048;                 //alk s mol/m^3
        y[4] = 3.052;                //alk deep mol/m^3
        y[5] = 0.367;                //po4 s mmol/m^3
        y[6] = 1.78;                //po4 deep mmol/m^3
        y[7] = -5.68;               //d13C atm
        y[8] = 3.11;                //d13C surf
        y[9] = 1.93;               //d13C deep
        y[10] = 0.231;               //bio C
        y[11] = -27.79;              //bio d13C
    */
    initCO2 = y[0] = 3.38;                 //atm pCO2
    y[1] = 3.08;               //sigC s mol/m^3
    y[2] = 3.23;                //sig C deep mol/m^3
    y[3] = 3.19;                 //alk s mol/m^3
    y[4] = 3.20;                //alk deep mol/m^3
    y[5] = 0.325;                //po4 s mmol/m^3
    y[6] = 1.576;                //po4 deep mmol/m^3
    y[7] = -4.7;               //d13C atm
    y[8] = 4.0;                //d13C surf
    y[9] = 3.1;               //d13C deep
    initbio = y[10] = 0.227;               //bio C
    y[11] = -28.1;              //bio d13C

    initsco3 = 0.200;

    //initialize sed. reservoirs
    carbeq(swtemp, y[1], y[3], sco3, spco2);
    fractionation();
    char sedfile[10] = "tmp00.txt";
    for (int counter2 = 0; sfdepths - counter2; counter2++)
    {
        if (counter2)
            sedfile[4] = '1';
        fin.open(sedfile);
        for (int counter = 0; sfcells - counter; counter++)
            fin >> y[nvars + 45 * counter2 + counter];
        for (int counter = 0; sfcells - counter; counter++)
            fin >> y[nvars + sfcells + 45 * counter2 + counter];
        fin.close();
        for (int counter = 0; sfcells - counter; counter++)
            y[nvars + sfcells * 2 + 45 * counter2 + counter] = y[8] - fraccal;
    }

    //calculate basin geometry
    double ddepth = 6000 / sfcells;
    double depth, uarea, barea, area;
    for (int i = 0; sfcells - i; i++) {
        sfres[i][0] = depth = double(i) * ddepth;
        sfres[i][1] = depth + ddepth;
        if (sfres[i][1] < 1000)                    //fractional areas
        {
            uarea = 0.15 * pow((depth) / 1000, 0.5);
            barea = 0.15 * pow((sfres[i][1]) / 1000, 0.5);
        }
        else if (depth < 1000)
        {
            uarea = 0.15 * std::pow((depth) / 1000, 0.5);
            barea = 0.85 * std::pow((sfres[i][1] - 1000) / 5000, 2) + 0.15;
        }
        else
        {
            uarea = 0.85 * std::pow((depth - 1000) / 5000, 2) + 0.15;
            barea = 0.85 * std::pow((sfres[i][1] - 1000) / 5000, 2) + 0.15;
        }
        sfres[i][2] = area = barea - uarea;       //sed res fractional area
        sfres[i][3] = area * 3.6e14 * sfthick;        //sed reservoir volumes m^3
    }

    hydro = 1.0;
    biology(y[0], y[10]);

    return;
}
//---------------------------------------------------------------------------
void equations(const double time, Vec_I_DP& y, Vec_IO_DP& dydx)
{
    update(y, time);
    burial(y, dydx);
    /*                                         //for sed spinup
       for(int i = 0; 10 - i; i++) dydx[i] = 0.0;
    */
    //atm C                                           //for simulations
    dydx[0] = (spco2 - y[0]) / etimesa * atmco2;
    dydx[0] += volcano + inject;
    dydx[0] -= assim - resp;
    dydx[0] /= atmco2;

    //surface Sigma C
    dydx[1] = (y[0] - spco2) / etimesa * atmco2;
    dydx[1] -= (1.0 + ioratio) * exprod;
    dydx[1] += exprod * shelfc * (1.0 - orgburial);
    dydx[1] += (y[2] - y[1]) * mixing;
    dydx[1] += sigcriv + rivorg;
    dydx[1] += shelfdiss;
    dydx[1] += doc + cpoc + (1 - fpocburial) * fpoc;
    dydx[1] /= vsurf;

    //deep Sigma C
    dydx[2] = (1.0 - orgburial) * exprod * (1 - shelfc);
    dydx[2] += caldiss;
    dydx[2] -= (y[2] - y[1]) * mixing;
    dydx[2] /= vdeep;

    //surface Alk
    dydx[3] = (y[4] - y[3]) * mixing;
    dydx[3] -= (2 * ioratio - 0.15) * exprod;
    dydx[3] -= 0.15 * exprod * shelfc * (1.0 - orgburial);
    dydx[3] += alkriv;
    dydx[3] += shelfdiss * 2.0;
    dydx[3] -= 0.15 * (doc + cpoc + (1 - fpocburial) * fpoc);
    dydx[3] /= vsurf;

    //deep Alk
    dydx[4] = 2.0 * caldiss;
    dydx[4] -= 0.15 * exprod * (1.0 - orgburial) * (1.0 - shelfc);
    dydx[4] -= (y[4] - y[3]) * mixing;
    dydx[4] /= vdeep;

    //surface PO4
    dydx[5] = (y[6] - y[5]) * mixing;
    dydx[5] -= exprod / cpratio;
    dydx[5] += exprod * shelfc * (1.0 - orgburial) / cpratio;
    dydx[5] += po4riv;
    dydx[5] /= vsurf;

    //deep PO4
    dydx[6] = (y[5] - y[6]) * mixing;
    dydx[6] += exprod * (1.0 - shelfc) * (1 - orgburial) / cpratio;
    dydx[6] /= vdeep;

    //atm d13C
    dydx[7] = spco2 * (y[8] - fracsa) - y[0] * y[7];
    dydx[7] /= etimesa;
    dydx[7] += volcano / atmco2 * c13volcano;
    dydx[7] += inject / atmco2 * c13inject;
    dydx[7] -= assim / atmco2 * (y[7] + biofrac);
    dydx[7] += resp / atmco2 * y[11];
    dydx[7] -= y[7] * dydx[0];
    dydx[7] /= y[0];

    //surface d13C
    dydx[8] = (y[0] * y[7] - spco2 * (y[8] - fracsa)) / etimesa * atmco2;
    dydx[8] += (y[2] * y[9] - y[1] * y[8]) * mixing;
    dydx[8] -= ((y[8] - fracorg) + ioratio * (y[8] - fraccal)) * exprod;
    dydx[8] += (y[8] - fracorg) * exprod * (1.0 - orgburial) * shelfc;
    dydx[8] += (rivorg * c13rivorg + rivcal * c13rivcal);
    dydx[8] += (doc + cpoc + (1 - fpocburial) * fpoc) * y[11];
    dydx[8] += shelfdiss * shelfdiss13C;
    dydx[8] /= vsurf;
    dydx[8] -= y[8] * dydx[1];
    dydx[8] /= y[1];

    //deep d13C
    dydx[9] = (y[1] * y[8] - y[2] * y[9]) * mixing;
    dydx[9] += (y[8] - fracorg) * exprod * (1 - orgburial) * (1.0 - shelfc);
    dydx[9] += caldiss * c13diss;
    dydx[9] /= vdeep;
    dydx[9] -= y[9] * dydx[2];
    dydx[9] /= y[2];

    //Bio C
    dydx[10] = assim;
    dydx[10] -= resp;
    dydx[10] -= eros;
    double dummy = dydx[10];
    //Bio d13C
    dydx[11] = assim * (y[7] + biofrac);
    dydx[11] -= (resp + eros) * y[11];
    dydx[11] -= y[11] * dydx[10];
    dydx[11] /= y[10];

    return;
}
//---------------------------------------------------------------------------
void update(Vec_I_DP& y, const double time) //update dependant variables
{
    if (inj) {
        injection(time);
        warm(time);
    }

    carbeq(swtemp, y[1], y[3], sco3, spco2);
    carbeq(dwtemp, y[2], y[4], dco3, dpco2);
    lysocline();
    production(y[5]);
    etimesa = 50;
    if (fb_weath)
    {
        rivcal = cwz * y[0];
        rivsil = swz * pow(y[0], 0.3);
        po4riv = pwz * pow(y[0], 0.3);
        rivorg = owz * pow(y[0], 0.3);
    }
    else
    {
        rivcal = cwz * initCO2;
        rivsil = swz * pow(initCO2, 0.3);
        po4riv = pwz * pow(initCO2, 0.3);
        rivorg = owz * pow(initCO2, 0.3);
    }
    alkriv = (2 * (rivcal + rivsil));
    sigcriv = rivcal;
    fractionation();
    if (fb_bio) biology(y[0], y[10]);
    else biology(initCO2, y[10]);

    return;
}
//---------------------------------------------------------------------------
void injection(const double time)
{
    double duration = 10000.0;
    double injrate = 0.3 / duration;

    if (ic == 0) {                       //first injection @ 20 ky
        if (time > 20000.0) {
            tw = ti = time;
            ic++;
        }
    }

    /*	if(ic > 0 && time - ti < 1000) inject = injrate;

        if(ic == 1 || ic == 2){  //second & third inj follow by 10 ky
            if(time - ti > 10000.0){
                    ti = time;
                    ic++;
            }
        }
    */
    //  Methane addition ramps linearly from 0 to injrate over 1000 ky

    if (ic > 0 && time - ti < 11000.0)
    {
        if (time - ti < 500)	inject = injrate * ((time - ti) / 500);
        else if (time - ti < 10500) inject = injrate;
        else if (time - ti < 11000 && ic != 3) inject = injrate *
            (11000 - (time - ti)) / 500;
    }

    else {
        inject = 0;
    }

}
//----------------------------------------------------------------------
void warm(const double time)
{
    //temperature as function of CO2 concentration and sensitivity
 /*	double sens = 5.0;
    swtemp = 290 + (pCO2 / 1.25 - 1.0) * sens;
    dwtemp = 281 + (pCO2 / 1.25 - 1.0) * sens;
                                                         */

                                                         //Perscribed temperature change for surface and deep ocean


    double a = 20000;
    double b = 30000;
    double c = 15000;

    double ab = a + b;
    double abc = a + b + c;

    if (time > 20000.0)
    {
        if (time - tw < a)
        {
            swtemp = 290 + 5 * (time - tw) / a;
            dwtemp = 281 + 5 * (time - tw) / a;
        }
        else if (time - tw > ab && time - tw < abc)
        {
            swtemp = 295 - 5 * (time - tw - ab) / c;
            dwtemp = 286 - 5 * (time - tw - ab) / c;
        }
        else if (time - tw > abc)
        {
            swtemp = 290;
            dwtemp = 281;
        }
    }

    //marine org carbon respiration as Q10 function
    if (fb_oc)
    {
        orgburial = 0.01 / pow(3.0, (swtemp - 290.0) / 10.0);
        fpocburial = 0.40 / pow(3.0, (swtemp - 290.0) / 10.0);
    }

    return;
}

//---------------------------------------------------------------------------
void carbeq(const double temp, const double sigc, const double alk,
    double& co3, double& pco2)
{
    double kcarb = 0.000575 + 0.000006 * (temp - 278);
    double kco2 = 0.035 + 0.0019 * (temp - 278);
    double hco3 = (sigc - sqrt(sigc * sigc - alk *
        (2 * sigc - alk) * (1 - 4 * kcarb))) /
        (1 - 4 * kcarb);
    co3 = (alk - hco3) / 2;
    pco2 = kco2 * hco3 * hco3 / co3;
    return;
}
//---------------------------------------------------------------------------
void lysocline()
{
    lyso = log(dco3 / 0.05779) / 0.16 + 4.0;
    if (lyso < 0.0) lyso = 0.0;
    if (lyso > 6.0) lyso = 6.0;
    return;
}
//---------------------------------------------------------------------------
void production(const double spo4)
{
    exprod = exratio * cpratio * spo4 * vsurf * prodfac;
    return;
}
//---------------------------------------------------------------------------
void fractionation()
{
    fracsa = 9483 / swtemp - 23.89;
    fraccal = 13 * (sco3 - initsco3);
    return;
}
//---------------------------------------------------------------------------
void biology(const double pco2, const double bioC)
{
    assim = 0.01 * pow(min(pco2, 2.5), 0.4);
    if (fb_bio) resp = bioC * 0.0634 * pow(3.0, (swtemp - 290.0) / 10.0);      //resp = 0.009967 * pow(pco2, 0.4)
    else resp = bioC * 0.0634;
    fpoc = 8.5e-6 * hydro;
    cpoc = 6.9e-6 * pow(bioC / 0.1537, 1.1);
    doc = 1.76e-5 * pow(bioC / 0.1537, 1.1);
    eros = fpoc + cpoc + doc;
    biofrac = -19 * pow(pco2, 0.17);

}
//---------------------------------------------------------------------------
void burial(Vec_I_DP& y, Vec_IO_DP& dydx)
{
    int ci, si, c13i, uci, usi, uc13i;
    double dissolution;              //carbonate dissolution
    double satu, satl, sat;
    double prezu, prezl, prez;
    double crain;               //carbonate flux to sediments
    double pores;
    double sigmaburial = 0.0;
    double sigmadiss = 0.0;
    double sigmac13diss = 0.0;
    double dissed, dissrain;
    double dc, ds;
    double sdc1, sdc2;
    double train = exprod * ioratio * 1.0e20;        //g calcite

    shelfc = sfres[0][2] * shelfratio;
    sdc1 = 1 - sfres[0][2];
    sdc2 = 1 - shelfc;
    for (int jrow = 0; sfcells - jrow; jrow++)
    {
        //indicies for sediment variables
        ci = nvars + jrow;
        si = nvars + sfcells + jrow;
        c13i = nvars + sfcells * 2 + jrow;

        //calculate rain accumulation and dissolution
        satu = dco3 / (0.05779 * exp(0.16 * (sfres[jrow][0] * 0.001 - 4.0)));
        if (satu > 1.0) satu = 1.0;
        satl = dco3 / (0.05779 * exp(0.16 * (sfres[jrow][1] * 0.001 - 4.0)));
        if (satl > 1.0) satl = 1.0;
        sat = (satu + satl) / 2.0;
        prezu = particle + (1.0 - particle)
            * exp(-(sfres[jrow][0] * 0.001 - lyso) / 0.25);
        if (prezu > 1.0) prezu = 1.0;
        prezl = particle + (1.0 - particle)
            * exp(-(sfres[jrow][1] * 0.001 - lyso) / 0.25);
        if (prezl > 1.0) prezl = 1.0;
        prez = (prezl + prezu) / 2.0;
        if (jrow == 0) {
            crain = shelfc * train * prez;
            dissrain = shelfc * train * (1.0 - prez);
        }
        else {
            crain = sfres[jrow][2] * train * prez * sdc2 / sdc1;
            dissrain = sfres[jrow][2] * train * (1.0 - prez) * sdc2 / sdc1;
        }
        //calculate sediment dissolution
        pores = porosity(y[ci], y[si], 0);
        if (1.0 - sat) {
            dissed = 360.0 * y[ci] * pow(1.0 - sat, 4.5) / (y[ci] + y[si]);
            dissed *= y[ci] * 0.1 * pores * pores;
        }
        else dissed = 0.0;
        if (jrow == 0) dissed += shelfc * train / ioratio * oxyresp;
        else dissed += sfres[jrow][2] * train / ioratio * oxyresp;
        if (dissed > y[ci]) dissed = y[ci];

        //sum dissolution and burial
        dissolution = dissed + dissrain;
        sigmaburial += (crain - dissed);
        if (jrow == 0) shelfdiss = dissolution;
        else sigmadiss += dissolution;
        if (jrow == 0) shelfdiss13C = dissed * y[c13i] + (y[8] - fraccal) * dissrain;
        else sigmac13diss += dissed * y[c13i] + (y[8] - fraccal) * dissrain;

        //initial derivatives for surface sediments
        dydx[si] = sfres[jrow][3] / 0.1 * 1.0;
        dydx[ci] = crain - dissed;
        dydx[c13i] = (y[8] - fraccal) * crain - y[c13i] * dissed;

        //calculate transfer due to inflation or deflation
        for (int jcol = 0; sfdepths - jcol; jcol++)
        {
            uci = ci, usi = si, uc13i = c13i;
            ci += 3 * sfcells;
            si += 3 * sfcells;
            c13i += 3 * sfcells;
            pores = porosity(y[uci], y[usi], jcol);

            //transfer between layers 1-5
            if (sfdepths - jcol - 1) {
                ds = Ds(sfres[jrow][3], y[uci], y[ci], y[usi], y[si], pores, dc, jcol);
                dydx[uci] += dc;
                dydx[usi] += ds;
                dydx[ci] = -dc;
                dydx[si] = -ds;
                if (ds >= 0.0) {
                    dydx[uc13i] += y[c13i] * dc;
                    dydx[uc13i] -= y[uc13i] * dydx[uci];
                    dydx[uc13i] /= y[uci];
                    dydx[c13i] = y[c13i] * -dc;
                }
                else {
                    dydx[uc13i] += y[uc13i] * dc;
                    dydx[uc13i] -= y[uc13i] * dydx[uci];
                    dydx[uc13i] /= y[uci];
                    dydx[c13i] = y[uc13i] * -dc;
                }
            }

            //layer 5 volume conservation
            else {
                ds = Ds(sfres[jrow][3], y[uci], y[uci], y[usi], y[usi], pores, dc, jcol);
                dydx[uci] += dc;
                dydx[usi] += ds;
                dydx[uc13i] += y[uc13i] * dc - y[uc13i] * dydx[uci];
                dydx[uc13i] /= y[uci];
            }
        }
    }
    calburial = sigmaburial / 1.0e20;     //10^18 mol
    caldiss = sigmadiss / 1.0e20;
    shelfdiss13C /= shelfdiss;
    shelfdiss /= 1.0e20;
    c13diss = sigmac13diss / sigmadiss;
    return;
}
//---------------------------------------------------------------------------
double porosity(double c, double s, int ind)
{
    double pores;
    double pctC = (c) / (c + s);
    double pmax = 1 - (0.483 + 0.45 * pctC) / 2.5;
    double alpha = 0.25 * pctC + 3.0 * (1 - pctC);
    pores = pmax + (1 - pmax) * exp(-((ind * sfthick + sfthick / 2) * 100 / alpha));
    return pores;
}
//--------------------------------------------------------------------------
double Ds(double vn, double c, double ct, double s, double st, double pores, double& dc, int ind)
{
    const double vc = 2720000.0;         //density, g/m^3
    const double vs = 2500000.0;

    double ds = vn * (1.0 - pores) - (c / vc + s / vs);
    if (ds > 0.0) {                             //defation
        double tpores;
        tpores = porosity(ct, st, ind);

        dc = ds * (ct / vc) / ((ct / vc + st / vs) / (1 - tpores));
        ds *= (st / vs) / ((ct / vc + st / vs) / (1 - tpores));
    }
    else {                                     //inflation
        dc = ds * (c / vc) / ((c / vc + s / vs) / (1 - pores));
        ds *= (s / vs) / ((c / vc + s / vs) / (1 - pores));
    }
    ds *= vs;
    dc *= vc;

    return ds;
}



