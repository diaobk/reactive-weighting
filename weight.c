#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double linear_interp(double, double, double, double, double);

double nonreactive_radial_force(double, double, double, double, double, double);

int main()
{
    //pseudocode
    //read in the tabulated force from the nonreactive model
    //read in the eigenavalues of the frame
    //read in forces and positions of each frame of the mapped reactive model
    //calculate the relavent distance/varibable for group of beads. Look up value in table and subtract away weighted value
    //	-MeOH beads and reactive chloride. Distance
    //	-MeOH beads and reactive carbon. Distance
    //	-MeOH beads and bonded chlorine. Distance
    //	-MeOH beads and propyl bead. Distance
    //	-cesium and reactive chloride. Distance
    //	-cesium and reactive carbon. Distance
    //	-cesium and bonded chlorine. Distance
    //	-cesium and propyl. Distance
    //	-reactive chloride and reactive carbon. Distance
    //	-reactive chloride and bonded chlorine. Distance
    //	-reactive chloride and propyl bead. Distance
    // 	-reactive carbon and bonded chloride. Distance
    //	-reactive carbon and propyl bead. Distance
    // 	-bonded chloride and propyl bead. Distance
    //	-reactive chloride, reactive carbon, and propyl bead. Angle
    //	-reactive carbon, propyl bead, and bonded chlorine. Angle
    //divide the remaining force by twice the product of the eigenvalues

    FILE *ifp;
    FILE *ifp0;
    FILE *ifp1;
    FILE *ifp2;
    FILE *ifp3;
    FILE *ifp4;
    FILE *ifp5;
    FILE *ifp6;
    FILE *ifp7;
    FILE *ifp8;
    FILE *ifp9;
    FILE *ifp10;
    FILE *ifp11;
    FILE *ifp12;
    FILE *ifp13;
    FILE *ifp14;
    FILE *ofp;

    ifp = fopen("MeOH_sample1_nobias.cg.dat", "r");
    ifp0 = fopen("coeff.dat", "r");
    ifp1 = fopen("fmecMeOH_Clion.table", "r");
    ifp2 = fopen("fmecMeOH_C1.table", "r");
    ifp3 = fopen("fmecMeOH_ClC.table", "r");
    ifp4 = fopen("fmecMeOH_C2.table", "r");
    ifp5 = fopen("fmecClion_Csion.table", "r");
    ifp6 = fopen("fmecCsion_C1.table", "r");
    ifp7 = fopen("fmecCsion_ClC.table", "r");
    ifp8 = fopen("fmecCsion_C2.table", "r");
    ifp9 = fopen("fmecClion_C1.table", "r");
    ifp10 = fopen("fmecClion_ClC.table", "r");
    ifp11 = fopen("fmecClion_C2.table", "r");
    ifp12 = fopen("fmecC1_ClC_bon.table", "r");
    ifp13 = fopen("fmecC1_C2_bon.table", "r");
    ifp14 = fopen("fmecMeOH_MeOH.table", "r");
    ofp = fopen("cg.MeOH_sample_remainder.dat", "w");

    const unsigned n_total_beads = 1005;
    const unsigned n_methanol_beads = 1000;
    const unsigned reactive_cl_id = 1000;
    const unsigned positive_ion_id = 1001;
    const unsigned reactive_c_id = 1002;
    const unsigned bonded_cl_id = 1003;
    const unsigned propyl_id = 1004;

    int mol[n_total_beads];
    int type[n_total_beads];
    double mass[n_total_beads];
    double q[n_total_beads];

    double x[n_total_beads];
    double y[n_total_beads];
    double z[n_total_beads];

    double fx[n_total_beads];
    double fy[n_total_beads];
    double fz[n_total_beads];

    double c[2];

    // Nonbonded force tables
    double f[12][201];
    const double nonbond_table_binwidth = 0.05;
    
    // Bonded force tables
    double fbond[2][408];
    const double bond_table_binwidth = 0.01;

    double disx;
    double disy;
    double disz;

    double disx2;
    double disy2;
    double disz2;

    double r;
    double r2;
    double s;
    double theta;

    int nonbond_table_bin;
    int bond_table_bin;

    double cgforce0;
    double cgforce1;

    double a;
    double a11;
    double a12;
    double a22;
    double dot;

    double nbox;
    double box;

    double sum;

    int t;
    int i;
    int j;
    int k;
    int id;
    int ts;

    char buffer[100];

    // read in cg force files
    //MeOH Clion
    for (i = 0; i < 201; i++) {
        fscanf(ifp1, "%lf %lf\n", &r, &f[0][i]);
    }
    //MeOH C1
    for (i = 0; i < 201; i++) {
        fscanf(ifp2, "%lf %lf\n", &r, &f[1][i]);
    }
    //MeOH ClC
    for (i = 0; i < 201; i++) {
        fscanf(ifp3, "%lf %lf\n", &r, &f[2][i]);
    }
    //MeOH C2
    for (i = 0; i < 201; i++) {
        fscanf(ifp4, "%lf %lf\n", &r, &f[3][i]);
    }
    //Csion Clion
    for (i = 0; i < 201; i++) {
        fscanf(ifp5, "%lf %lf\n", &r, &f[4][i]);
    }
    //Csion C1
    for (i = 0; i < 201; i++) {
        fscanf(ifp6, "%lf %lf\n", &r, &f[5][i]);
    }
    //Csion ClC
    for (i = 0; i < 201; i++) {
        fscanf(ifp7, "%lf %lf\n", &r, &f[6][i]);
    }
    //Csion C2
    for (i = 0; i < 201; i++) {
        fscanf(ifp8, "%lf %lf\n", &r, &f[7][i]);
    }
    //Clion C1
    for (i = 0; i < 201; i++) {
        fscanf(ifp9, "%lf %lf\n", &r, &f[8][i]);
    }
    //Clion ClC
    for (i = 0; i < 201; i++) {
        fscanf(ifp10, "%lf %lf\n", &r, &f[9][i]);
    }
    //Clion C2
    for (i = 0; i < 201; i++) {
        fscanf(ifp11, "%lf %lf\n", &r, &f[10][i]);
    }
    //C1 ClC
    for (i = 0; i < 408; i++) {
        fscanf(ifp12, "%lf %lf\n", &r, &fbond[0][i]);
    }
    //C1 C2
    for (i = 0; i < 408; i++) {
        fscanf(ifp13, "%lf %lf\n", &r, &fbond[1][i]);
    }
    //MeOH MeOH
    for (i = 0; i < 201; i++) {
        fscanf(ifp14, "%lf %lf\n", &r, &f[11][i]);
    }

    //t is timestep counter

    for (t = 0; t < 1; t++) {

        //read in line of coefficient file

        fscanf(ifp0, "%lf %lf\n", &c[0], &c[1]);

        printf("%lf %lf\n", c[0], c[1]);

        //read in header

        printf("%d\n", t);

        fgets(buffer, 100, ifp);
        fprintf(ofp, "ITEM: TIMESTEP\n");

        fscanf(ifp, "%d\n", &ts);
        fprintf(ofp, "%d\n", ts);

        fgets(buffer, 100, ifp);
        fprintf(ofp, "ITEM: NUMBER OF ATOMS\n");

        fgets(buffer, 100, ifp);
        fprintf(ofp, "1689\n");

        fgets(buffer, 100, ifp);
        fprintf(ofp, "ITEM: BOX BOUNDS pp pp pp\n");

        fscanf(ifp, "%lf %lf\n", &nbox, &box);
        fprintf(ofp, "%lf %lf\n", nbox, box);

        fgets(buffer, 100, ifp);
        fprintf(ofp, "%lf %lf\n", nbox, box);

        fgets(buffer, 100, ifp);
        fprintf(ofp, "%lf %lf\n", nbox, box);

        fgets(buffer, 100, ifp);
        fprintf(ofp, "ITEM: ATOMS id mol type x y z fx fy fz \n");
        //printf("%s\n",buffer);
        printf("%f\n", box - nbox);
        //printf("%f\n",nbox);

        for (i = 0; i < n_total_beads; i++) {

            //read in coordinates

            fscanf(ifp, "%d %d %d %lf %lf %lf %lf %lf %lf\n", &id, &mol[i], &type[i], &x[i], &y[i], &z[i], &fx[i], &fy[i], &fz[i]);
        }
        /*
        		for(i=0;i<999;i++){
        			for(j=i+1;j<n_methanol_beads;j++){

        				disx = x[i]-x[j];
        				disy = y[i]-y[j];
        				disz = z[i]-z[j];

        				if(disx > (box-nbox)/2.0){
        					disx = disx - (box-nbox);
        				}
        				if(disx < -(box-nbox)/2.0){
        					disx = disx + (box-nbox);
        				}
        				if(disy > (box-nbox)/2.0){
        					disy = disy - (box-nbox);
        				}
        				if(disy < -(box-nbox)/2.0){
        					disy = disy + (box-nbox);
        				}
        				if(disz > (box-nbox)/2.0){
        					disz = disz - (box-nbox);
        				}
        				if(disz < -(box-nbox)/2.0){
        					disz = disz + (box-nbox);
        				}

        				r = sqrt(disx*disx + disy*disy + disz*disz);


        				if( r < 10.0){
        				nonbond_table_bin = r/nonbond_table_binwidth;

        //				printf("%lf %d\n",r,nonbond_table_bin);

        				cgforce0 = linear_interp(nonbond_table_bin*nonbond_table_binwidth,(nonbond_table_binwidth*nonbond_table_bin+nonbond_table_binwidth),f[11][nonbond_table_bin],f[11][nonbond_table_bin+1],r);

        //				printf("%lf\n",cgforce0);

        				fx[j] = fx[j] + (cgforce0*disx/r);
        				fy[j] = fy[j] + (cgforce0*disy/r);
        				fz[j] = fz[j] + (cgforce0*disz/r);

        				fx[i] = fx[i] - (cgforce0*disx/r);
        				fy[i] = fy[i] - (cgforce0*disy/r);
        				fz[i] = fz[i] - (cgforce0*disz/r);
        				}
        			}
        		}
        */
        for (i = 0; i < n_methanol_beads; i++) {

            //calculate distance between methanol and reactive chloride.

            disx = x[i] - x[reactive_cl_id];
            disy = y[i] - y[reactive_cl_id];
            disz = z[i] - z[reactive_cl_id];

            //periodic boundary conditions

            if (disx > (box - nbox) / 2.0) {
                disx = disx - (box - nbox);
            }
            if (disx < -(box - nbox) / 2.0) {
                disx = disx + (box - nbox);
            }
            if (disy > (box - nbox) / 2.0) {
                disy = disy - (box - nbox);
            }
            if (disy < -(box - nbox) / 2.0) {
                disy = disy + (box - nbox);
            }
            if (disz > (box - nbox) / 2.0) {
                disz = disz - (box - nbox);
            }
            if (disz < -(box - nbox) / 2.0) {
                disz = disz + (box - nbox);
            }

            r = sqrt(disx * disx + disy * disy + disz * disz);

//			fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r);
//			fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r);
//			fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r);

            //total reactive chloride force = total reactive choride force - coefficient * coeffcient * table value for reactive chloride and methanol - coeffcient *coeffcient * table value for bonded chloride and methanol
            if (r < 10.0) {
                nonbond_table_bin = r / nonbond_table_binwidth;

//			printf("%lf %d\n",r,nonbond_table_bin);

                cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_binwidth * nonbond_table_bin + nonbond_table_binwidth), f[0][nonbond_table_bin], f[0][nonbond_table_bin + 1], r);
                cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_binwidth * nonbond_table_bin + nonbond_table_binwidth), f[2][nonbond_table_bin], f[2][nonbond_table_bin + 1], r);

//			printf("%lf %lf\n",cgforce0,cgforce1);


                fx[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
                fy[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
                fz[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
            }

///////////////////////////////////////////

            //methanol and reactive carbon bead

            disx = x[i] - x[reactive_c_id];
            disy = y[i] - y[reactive_c_id];
            disz = z[i] - z[reactive_c_id];

            if (disx > (box - nbox) / 2.0) {
                disx = disx - (box - nbox);
            }
            if (disx < -(box - nbox) / 2.0) {
                disx = disx + (box - nbox);
            }
            if (disy > (box - nbox) / 2.0) {
                disy = disy - (box - nbox);
            }
            if (disy < -(box - nbox) / 2.0) {
                disy = disy + (box - nbox);
            }
            if (disz > (box - nbox) / 2.0) {
                disz = disz - (box - nbox);
            }
            if (disz < -(box - nbox) / 2.0) {
                disz = disz + (box - nbox);
            }

            r = sqrt(disx * disx + disy * disy + disz * disz);

//			fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r);
//			fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r);
//			fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r);

            //total reactive carbon force = total reactive carbon force - coefficient * coeffcient * table value for reactive carbon and methanol - coefficient * coeffcient * table value for reactive carbon and methanol
            if (r < 10.0) {
                nonbond_table_bin = r / nonbond_table_binwidth;

                cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[1][nonbond_table_bin], f[1][nonbond_table_bin + 1], r);
                cgforce1 = cgforce0;

                fx[reactive_c_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
                fy[reactive_c_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
                fz[reactive_c_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
            }
///////////////////////////////////////////

            //calculate distance between methanol and bonded chloride

            disx = x[i] - x[bonded_cl_id];
            disy = y[i] - y[bonded_cl_id];
            disz = z[i] - z[bonded_cl_id];

            if (disx > (box - nbox) / 2.0) {
                disx = disx - (box - nbox);
            }
            if (disx < -(box - nbox) / 2.0) {
                disx = disx + (box - nbox);
            }
            if (disy > (box - nbox) / 2.0) {
                disy = disy - (box - nbox);
            }
            if (disy < -(box - nbox) / 2.0) {
                disy = disy + (box - nbox);
            }
            if (disz > (box - nbox) / 2.0) {
                disz = disz - (box - nbox);
            }
            if (disz < -(box - nbox) / 2.0) {
                disz = disz + (box - nbox);
            }

            r = sqrt(disx * disx + disy * disy + disz * disz);

//			fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//			fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//			fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r)

            //total bonded chlorine force = total bonded chlorine force - coefficient * coeffcient * table value for bonded chlorine and methanol - coeffcient *coeffcient * table value for reactive chlorine  and methanol
            if (r < 10.0) {
                nonbond_table_bin = r / nonbond_table_binwidth;

                cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[2][nonbond_table_bin], f[2][nonbond_table_bin + 1], r);
                cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[0][nonbond_table_bin], f[0][nonbond_table_bin + 1], r);

                fx[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
                fy[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
                fz[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
            }
/////////////////////////////////////////////

            //methanol and propyl bead

            disx = x[i] - x[propyl_id];
            disy = y[i] - y[propyl_id];
            disz = z[i] - z[propyl_id];

            if (disx > (box - nbox) / 2.0) {
                disx = disx - (box - nbox);
            }
            if (disx < -(box - nbox) / 2.0) {
                disx = disx + (box - nbox);
            }
            if (disy > (box - nbox) / 2.0) {
                disy = disy - (box - nbox);
            }
            if (disy < -(box - nbox) / 2.0) {
                disy = disy + (box - nbox);
            }
            if (disz > (box - nbox) / 2.0) {
                disz = disz - (box - nbox);
            }
            if (disz < -(box - nbox) / 2.0) {
                disz = disz + (box - nbox);
            }

            r = sqrt(disx * disx + disy * disy + disz * disz);

//			fx[i] -= c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//			fy[i] -= c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//			fz[i] -= c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r)

            //total propyl force = total propyl  force - coefficient * coeffcient * table value for propyl and methanol - coefficient * coeffcient * table value for propyl and methanol
            if (r < 10.0) {
                nonbond_table_bin = r / nonbond_table_binwidth;

                cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[3][nonbond_table_bin], f[3][nonbond_table_bin + 1], r);
                cgforce1 = cgforce0;

                fx[propyl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
                fy[propyl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
                fz[propyl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
            }
        }

///////////////////////////////////////////

        //positive ion and reactive chloride

        disx = x[positive_ion_id] - x[reactive_cl_id];
        disy = y[positive_ion_id] - y[reactive_cl_id];
        disz = z[positive_ion_id] - z[reactive_cl_id];

        //periodic boundary conditions

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

//		fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r)

        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[4][nonbond_table_bin], f[4][nonbond_table_bin + 1], r);
            cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[6][nonbond_table_bin], f[6][nonbond_table_bin + 1], r);

            fx[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //positive ion and reactive carbon bead

        disx = x[positive_ion_id] - x[reactive_c_id];
        disy = y[positive_ion_id] - y[reactive_c_id];
        disz = z[positive_ion_id] - z[reactive_c_id];

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }
        r = sqrt(disx * disx + disy * disy + disz * disz);

//		fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r)

        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_binwidth * nonbond_table_bin + nonbond_table_binwidth), f[5][nonbond_table_bin], f[5][nonbond_table_bin + 1], r);
            cgforce1 = cgforce0;

            fx[reactive_c_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_c_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_c_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //positive ion and bonded chlorine

        disx = x[positive_ion_id] - x[bonded_cl_id];
        disy = y[positive_ion_id] - y[bonded_cl_id];
        disz = z[positive_ion_id] - z[bonded_cl_id];

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

//		fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*di[0]sz/
        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[6][nonbond_table_bin], f[6][nonbond_table_bin + 1], r);
            cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[4][nonbond_table_bin], f[4][nonbond_table_bin + 1], r);

            fx[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //positive ion and propyl bead

        disx = x[positive_ion_id] - x[propyl_id];
        disy = y[positive_ion_id] - y[propyl_id];
        disz = z[positive_ion_id] - z[propyl_id];

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

//		fx[i] = fx[i] - c[t][0]*(f[0][r]*disx/r)-c[t][1](f[0][r]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[0][r]*disy/r)-c[t][1](f[0][r]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[0][r]*disz/r)-c[t][1](f[0][r]*disz/r)

        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[7][nonbond_table_bin], f[7][nonbond_table_bin + 1], r);
            cgforce1 = cgforce0;


            fx[propyl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[propyl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[propyl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //reactive chloride with reactive carbon bead

        disx = x[reactive_cl_id] - x[reactive_c_id];
        disy = y[reactive_cl_id] - y[reactive_c_id];
        disz = z[reactive_cl_id] - z[reactive_c_id];

        //periodic boundary conditions

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;
            bond_table_bin = r / bond_table_binwidth;
            
            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[8][nonbond_table_bin], f[8][nonbond_table_bin + 1], r);
            cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, fbond[1][bond_table_bin], fbond[1][bond_table_bin + 1], r);

            fx[reactive_c_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_c_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_c_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);

            fx[reactive_cl_id] += - c[0] * c[0] * (cgforce0 * disx / r) - c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_cl_id] += - c[0] * c[0] * (cgforce0 * disy / r) - c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_cl_id] += - c[0] * c[0] * (cgforce0 * disz / r) - c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //reactive chloride with bonded chloride

        disx = x[reactive_cl_id] - x[bonded_cl_id];
        disy = y[reactive_cl_id] - y[bonded_cl_id];
        disz = z[reactive_cl_id] - z[bonded_cl_id];

        //periodic boundary conditions

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);


        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[9][nonbond_table_bin], f[9][nonbond_table_bin + 1], r);
            cgforce1 = cgforce0;

            fx[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);

            fx[reactive_cl_id] -= c[0] * c[0] * (cgforce0 * disx / r) - c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_cl_id] -= c[0] * c[0] * (cgforce0 * disy / r) - c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_cl_id] -= c[0] * c[0] * (cgforce0 * disz / r) - c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////
        /*
        		//reactive chloride with propyl bead

        		disx = x[reactive_cl_id]-x[propyl_id];
        		disy = y[reactive_cl_id]-y[propyl_id];
        		disz = z[reactive_cl_id]-z[propyl_id];

        		//periodic boundary conditions

        		if(disx > (box-nbox)/2.0){
        			disx = disx - (box-nbox);
        		}
        		if(disx < -(box-nbox)/2.0){
        			disx = disx + (box-nbox);
        		}
        		if(disy > (box-nbox)/2.0){
        			disy = disy - (box-nbox);
        		}
        		if(disy < -(box-nbox)/2.0){
        			disy = disy + (box-nbox);
        		}
        		if(disz > (box-nbox)/2.0){
        			disz = disz - (box-nbox);
        		}
        		if(disz < -(box-nbox)/2.0){
        			disz = disz + (box-nbox);
        		}

        		r = sqrt(disx*disx + disy*disy + disz*disz);


        		if( r < 10.0){
        		nonbond_table_bin = r/nonbond_table_binwidth;

        		cgforce0 = linear_interp(nonbond_table_bin,(nonbond_table_bin+nonbond_table_binwidth),f[9][nonbond_table_bin],f[9][nonbond_table_bin+1],r);

        		fx[reactive_cl_id] = fx[reactive_cl_id] + c[0]*c[0]*(cgforce0*disx/r);
        		fy[reactive_cl_id] = fy[reactive_cl_id] + c[0]*c[0]*(cgforce0*disy/r);
        		fz[reactive_cl_id] = fz[reactive_cl_id] + c[0]*c[0]*(cgforce0*disz/r);

        		fx[reactive_c_id] = fx[reactive_c_id] - c[0]*c[0]*(cgforce0*disx/r);
        		fy[reactive_c_id] = fy[reactive_c_id] - c[0]*c[0]*(cgforce0*disy/r);
        		fz[reactive_c_id] = fz[reactive_c_id] - c[0]*c[0]*(cgforce0*disz/r);
        		}
        */
///////////////////////////////////////////

        //reactive carbon and bonded chloride

        disx = x[reactive_c_id] - x[bonded_cl_id];
        disy = y[reactive_c_id] - y[bonded_cl_id];
        disz = z[reactive_c_id] - z[bonded_cl_id];

        //periodic boundary conditions

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

        if (r < 10.0) {
            nonbond_table_bin = r / nonbond_table_binwidth;
            bond_table_bin = r / bond_table_binwidth;

            cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, fbond[1][bond_table_bin], fbond[1][bond_table_bin + 1], r);
            cgforce1 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, f[8][nonbond_table_bin], f[8][nonbond_table_bin + 1], r);

            fx[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
            fy[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
            fz[bonded_cl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);

            fx[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disx / r) - c[1] * c[1] * (cgforce1 * disx / r);
            fy[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disy / r) - c[1] * c[1] * (cgforce1 * disy / r);
            fz[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disz / r) - c[1] * c[1] * (cgforce1 * disz / r);
        }
///////////////////////////////////////////

        //reactive carbon and propyl

        disx = x[reactive_c_id] - x[propyl_id];
        disy = y[reactive_c_id] - y[propyl_id];
        disz = z[reactive_c_id] - z[propyl_id];

        //periodic boundary conditions

        if (disx > (box - nbox) / 2.0) {
            disx = disx - (box - nbox);
        }
        if (disx < -(box - nbox) / 2.0) {
            disx = disx + (box - nbox);
        }
        if (disy > (box - nbox) / 2.0) {
            disy = disy - (box - nbox);
        }
        if (disy < -(box - nbox) / 2.0) {
            disy = disy + (box - nbox);
        }
        if (disz > (box - nbox) / 2.0) {
            disz = disz - (box - nbox);
        }
        if (disz < -(box - nbox) / 2.0) {
            disz = disz + (box - nbox);
        }

        r = sqrt(disx * disx + disy * disy + disz * disz);

        nonbond_table_bin = r / nonbond_table_binwidth;

        cgforce0 = linear_interp(nonbond_table_bin * nonbond_table_binwidth, (nonbond_table_bin + 1) * nonbond_table_binwidth, fbond[1][nonbond_table_bin], fbond[1][nonbond_table_bin + 1], r);
        cgforce1 = cgforce0;

        fx[propyl_id] += c[0] * c[0] * (cgforce0 * disx / r) + c[1] * c[1] * (cgforce1 * disx / r);
        fy[propyl_id] += c[0] * c[0] * (cgforce0 * disy / r) + c[1] * c[1] * (cgforce1 * disy / r);
        fz[propyl_id] += c[0] * c[0] * (cgforce0 * disz / r) + c[1] * c[1] * (cgforce1 * disz / r);

        fx[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disx / r) - c[1] * c[1] * (cgforce1 * disx / r);
        fy[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disy / r) - c[1] * c[1] * (cgforce1 * disy / r);
        fz[reactive_c_id] -= c[0] * c[0] * (cgforce0 * disz / r) - c[1] * c[1] * (cgforce1 * disz / r);

///////////////////////////////////////////
        /*
        		//bonded chloride and propyl

        		disx = x[bonded_cl_id]-x[propyl_id];
        		disy = y[bonded_cl_id]-y[propyl_id];
        		disz = z[bonded_cl_id]-z[propyl_id];

        		//periodic boundary conditions

        		if(disx > (box-nbox)/2.0){
        			disx = disx - (box-nbox);
        		}
        		if(disx < -(box-nbox)/2.0){
        			disx = disx + (box-nbox);
        		}
        		if(disy > (box-nbox)/2.0){
        			disy = disy - (box-nbox);
        		}
        		if(disy < -(box-nbox)/2.0){
        			disy = disy + (box-nbox);
        		}
        		if(disz > (box-nbox)/2.0){
        			disz = disz - (box-nbox);
        		}
        		if(disz < -(box-nbox)/2.0){
        			disz = disz + (box-nbox);
        		}

        		r = sqrt(disx*disx + disy*disy + disz*disz);

        		if( r < 10.0){

        		nonbond_table_bin=r/nonbond_table_binwidth;

        		cgforce1 = linear_interp(nonbond_table_bin,(nonbond_table_bin+nonbond_table_binwidth),fbond[1][nonbond_table_bin],fbond[1][nonbond_table_bin+1],r)


        		fx[bonded_cl_id] -= c[1]*c[1]*(cgforce1*disx/r);
        		fy[bonded_cl_id] -= c[1]*c[1]*(cgforce1*disy/r);
        		fz[bonded_cl_id] -= c[1]*c[1]*(cgforce1*disz/r);

        		fx[propyl_id] += c[1]*c[1]*(cgforce1*disx/r);
        		fy[propyl_id] += c[1]*c[1]*(cgforce1*disy/r);
        		fz[propyl_id] += c[1]*c[1]*(cgforce1*disz/r);
        		}
        */
///////////////////////////////////////////
        /*

        		//angle with reactive

        		disx = x[reactive_cl_id]-x[reactive_c_id];
        		disy = y[reactive_cl_id]-y[reactive_c_id];
        		disz = z[reactive_cl_id]-z[reactive_c_id];

        		disx2 = x[propyl_id]-x[reactive_c_id];
        		disy2 = y[propyl_id]-y[reactive_c_id];
        		disz2 = z[propyl_id]-z[reactive_c_id];

        		if(disx > (box-nbox)/2.0){
        			disx = disx - (box-nbox);
        		}
        		if(disx < -(box-nbox)/2.0){
        			disx = disx + (box-nbox);
        		}
        		if(disy > (box-nbox)/2.0){
        			disy = disy - (box-nbox);
        		}
        		if(disy < -(box-nbox)/2.0){
        			disy = disy + (box-nbox);
        		}
        		if(disz > (box-nbox)/2.0){
        			disz = disz - (box-nbox);
        		}
        		if(disz < -(box-nbox)/2.0){
        			disz = disz + (box-nbox);
        		}

        		if(disx2 > (box-nbox)/2.0){
        			disx2 = disx2 - (box-nbox);
        		}
        		if(disx2 < -(box-nbox)/2.0){
        			disx2 = disx2 + (box-nbox);
        		}
        		if(disy2 > (box-nbox)/2.0){
        			disy2 = disy2 - (box-nbox);
        		}
        		if(disy2 < -(box-nbox)/2.0){
        			disy2 = disy2 + (box-nbox);
        		}
        		if(disz2 > (box-nbox)/2.0){
        			disz2 = disz2 - (box-nbox);
        		}
        		if(disz2 < -(box-nbox)/2.0){
        			disz2 = disz2 + (box-nbox);
        		}

        		r = sqrt(disx*disx + disy*disy + disz*disz);

        		r2 = sqrt(disx2*disx2 + disy2*disy2 + disz2*disz2);

        		dot = delx1*delx2 + dely1*dely2 + delz1*delz2;
        		dot /= r1*r2;
        		s = sqrt(1.0 - dot*dot);

        		theta = acos(dot);

        		a = f[theta][type] * s;
        		a11 = a*dot / (r*r);
        		a12 = -a / (r*r2);
        		a22 = a*dot / (r2*r3);

        		fx[reactive_cl_id] -= c[0]*c[0]*(a11*delx + a12*delx2);
        		fy[reactive_cl_id] -= c[0]*c[0]*(a11*dely + a12*dely2);
        		fz[reactive_cl_id] -= c[0]*c[0]*(a11*delz + a12*delz2);

        		fx[reactive_c_id] += c[0]*c[0]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
        		fy[reactive_c_id] += c[0]*c[0]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
        		fz[reactive_c_id] += c[0]*c[0]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

        		fx[propyl_id] -= c[0]*c[0]*(a22*delx2 + a12*delx);
        		fy[propyl_id] -= c[0]*c[0]*(a22*dely2 + a12*dely);
        		fz[propyl_id] -= c[0]*c[0]*(a22*delz2 + a12*delz);

        ////////////////////////////////

        		//angle with bonded

        		disx = x[bonded_cl_id]-x[reactive_c_id];
        		disy = y[bonded_cl_id]-y[reactive_c_id];
        		disz = z[bonded_cl_id]-z[reactive_c_id];

        		disx2 = x[propyl_id]-x[reactive_c_id];
        		disy2 = y[propyl_id]-y[reactive_c_id];
        		disz2 = z[propyl_id]-z[reactive_c_id];

        		if(disx > (box-nbox)/2.0){
        			disx = disx - (box-nbox);
        		}
        		if(disx < -(box-nbox)/2.0){
        			disx = disx + (box-nbox);
        		}
        		if(disy > (box-nbox)/2.0){
        			disy = disy - (box-nbox);
        		}
        		if(disy < -(box-nbox)/2.0){
        			disy = disy + (box-nbox);
        		}
        		if(disz > (box-nbox)/2.0){
        			disz = disz - (box-nbox);
        		}
        		if(disz < -(box-nbox)/2.0){
        			disz = disz + (box-nbox);
        		}

        		if(disx2 > (box-nbox)/2.0){
        			disx2 = disx2 - (box-nbox);
        		}
        		if(disx2 < -(box-nbox)/2.0){
        			disx2 = disx2 + (box-nbox);
        		}
        		if(disy2 > (box-nbox)/2.0){
        			disy2 = disy2 - (box-nbox);
        		}
        		if(disy2 < -(box-nbox)/2.0){
        			disy2 = disy2 + (box-nbox);
        		}
        		if(disz2 > (box-nbox)/2.0){
        			disz2 = disz2 - (box-nbox);
        		}
        		if(disz2 < -(box-nbox)/2.0){
        			disz2 = disz2 + (box-nbox);
        		}

        		r = sqrt(disx*disx + disy*disy + disz*disz);

        		r2 = sqrt(disx2*disx2 + disy2*disy2 + disz2*disz2);

        		dot = delx1*delx2 + dely1*dely2 + delz1*delz2;
        		dot /= r1*r2;
        		s = sqrt(1.0 - dot*dot);

        		theta = acos(dot);

        		a = f[theta][type] * s;
        		a11 = a*dot / (r*r);
        		a12 = -a / (r*r2);
        		a22 = a*dot / (r2*r3);

        		fx[bonded_cl_id] = fx[bonded_cl_id] - c[1]*c[1]*(a11*delx + a12*delx2);
        		fy[bonded_cl_id] = fy[bonded_cl_id] - c[1]*c[1]*(a11*dely + a12*dely2);
        		fz[bonded_cl_id] = fz[bonded_cl_id] - c[1]*c[1]*(a11*delz + a12*delz2);

        		fx[reactive_c_id] = fx[reactive_c_id] + c[1]*c[1]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
        		fy[reactive_c_id] = fy[reactive_c_id] + c[1]*c[1]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
        		fz[reactive_c_id] = fz[reactive_c_id] + c[1]*c[1]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

        		fx[propyl_id] = fx[propyl_id] - c[1]*c[1]*(a22*delx2 + a12*delx);
        		fy[propyl_id] = fy[propyl_id] - c[1]*c[1]*(a22*dely2 + a12*dely);
        		fz[propyl_id] = fz[propyl_id] - c[1]*c[1]*(a22*delz2 + a12*delz);
        */

        // Weight forces by the product of the state weights.

        fx[reactive_cl_id] /= (2 * c[0] * c[1]);
        fy[reactive_cl_id] /= (2 * c[0] * c[1]);
        fz[reactive_cl_id] /= (2 * c[0] * c[1]);

        fx[reactive_c_id] /= (2 * c[0] * c[1]);
        fy[reactive_c_id] /= (2 * c[0] * c[1]);
        fz[reactive_c_id] /= (2 * c[0] * c[1]);

        fx[bonded_cl_id] /= (2 * c[0] * c[1]);
        fy[bonded_cl_id] /= (2 * c[0] * c[1]);
        fz[bonded_cl_id] /= (2 * c[0] * c[1]);

        fx[propyl_id] /= (2 * c[0] * c[1]);
        fy[propyl_id] /= (2 * c[0] * c[1]);
        fz[propyl_id] /= (2 * c[0] * c[1]);

        for (i = 0; i < n_total_beads; i++) {
            fprintf(ofp, "%d %d %d %6g %6g %6g %6g %6g %6g\n", i + 1, mol[i], type[i], x[i], y[i], z[i], fx[i], fy[i], fz[i]);
        }
    }

    return 0;
}

// Find the value of a function y at point x by linearly
// interpolating between two points x0 and x1 of a 
// function with values y0 and y1 at those points.
double linear_interp(double x0, double x1, double y0, double y1, double x)
{
    double slope = (y1 - y0) / (x1 - x0);
    return (y0 + (x - x0) * slope);
}

// General force weighting for the diagonal contribution
// of a radial force.
double nonreactive_radial_force(double c0, double c1, double r, double dis, double cgforce0, double cgforce1)
{
    return c0 * c0 * (cgforce0 * dis / r) + c1 * c1 * (cgforce1 * dis / r);
}
