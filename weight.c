#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double linear_interp(double,double,double,double,double);

double nonreactive_radial_force(double,double,double,double,double,double);

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

	ifp = fopen("MeOH_sample.dat","r");
	ifp2 = fopen("MeOH_sample.dat","r");
	ofp = fopen("cg.MeOH_sample.dat","w");

	int mol[1005];
	int type[1005];
	double mass[1005];
	double q[1005];

	double x[1005];
	double y[1005];
	double z[1005];

	double fx[1005];
	double fy[1005];
	double fz[1005];

	double c[2];

	double f[200][11];
	double fbond[408][2];

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

	int r0;
	int r1;

	double cgforce1;
	double cgforce2;

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
	for(i=0;i<200;i++){fscanf(ifp1,"%lf %lf\n",r,f[i][0]);}
	//MeOH C1
	for(i=0;i<200;i++){fscanf(ifp2,"%lf %lf\n",r,f[i][1]);}
	//MeOH ClC
	for(i=0;i<200;i++){fscanf(ifp3,"%lf %lf\n",r,f[i][2]);}
	//MeOH C2
	for(i=0;i<200;i++){fscanf(ifp4,"%lf %lf\n",r,f[i][3]);}
	//Csion Clion
	for(i=0;i<200;i++){fscanf(ifp5,"%lf %lf\n",r,f[i][4]);}
	//Csion C1
	for(i=0;i<200;i++){fscanf(ifp6,"%lf %lf\n",r,f[i][5]);}
	//Csion ClC
	for(i=0;i<200;i++){fscanf(ifp7,"%lf %lf\n",r,f[i][6]);}
	//Csion C2
	for(i=0;i<200;i++){fscanf(ifp8,"%lf %lf\n",r,f[i][7]);}
	//Clion C1
	for(i=0;i<200;i++){fscanf(ifp9,"%lf %lf\n",r,f[i][8]);}
	//Clion ClC
	for(i=0;i<200;i++){fscanf(ifp10,"%lf %lf\n",r,f[i][9]);}
	//Clion C2
	for(i=0;i<200;i++){fscanf(ifp11,"%lf %lf\n",r,f[i][10]);}
	//C1 ClC
	for(i=0;i<408;i++){fscanf(ifp12,"%lf %lf\n",r,fbond[i][0]);}
	//C1 C2
	for(i=0;i<408;i++){fscanf(ifp12,"%lf %lf\n",r,fbond[i][1]);}



	//t is timestep counter

	for(t=0;t<8000;t++){

		//read in line of coefficient file

		fscanf(ifp2,"%lf %lf\n",&c[0],&c[1]);

		//read in header

		printf("%d\n",t);

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: TIMESTEP\n",buffer);
		
		fscanf(ifp,"%d\n",&ts);
		fprintf(ofp,"%d\n",ts);

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: NUMBER OF ATOMS\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"1689\n");

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: BOX BOUNDS pp pp pp\n");

		fscanf(ifp,"%lf %lf\n",&nbox,&box);
		fprintf(ofp,"%lf %lf\n",nbox,box);

		fgets(buffer,100,ifp);
		fprintf(ofp,"%lf %lf\n",nbox,box);

		fgets(buffer,100,ifp);
		fprintf(ofp,"%lf %lf\n",nbox,box);

		fgets(buffer,100,ifp);
		fprintf(ofp,"ITEM: ATOMS id mol type x y z fx fy fz \n");
		//printf("%s\n",buffer);		
		printf("%f\n",box-nbox);
		//printf("%f\n",nbox);

		for(i=0;i<1005;i++){

			//read in coordinates

			fscanf(ifp,"%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&id,&mol[j],&type[j],&q[j],&mass[j],&x[j],&y[j],&z[j],&fx[j],&fy[j],&fz[j]);
		}

		for(i=0;i<1000;i++){

			//calculate distance between methanol and reactive chloride. 

			disx = x[i]-x[1000];
			disy = y[i]-y[1000];
			disz = z[i]-z[1000];

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

			

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r);
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r);
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r);

			//total reactive chloride force = total reactive choride force - coefficient * coeffcient * table value for reactive chloride and methanol - coeffcient *coeffcient * table value for bonded chloride and methanol
			if( r < 10.0){
			r0 = r/0.05;

			cgforce1 = linear_interp(r0,(r0+0.05),f[r0][0],f[r0+1][0],r);
			cgforce2 = linear_interp(r0,(r0+0.05),f[r0][2],f[r0+1][2],r);


			fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
			fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
			fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
			}

///////////////////////////////////////////

			//methanol and reactive carbon bead

			disx = x[i]-x[1002];
			disy = y[i]-y[1002];
			disz = z[i]-z[1002];

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

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r);
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r);
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r);

			//total reactive carbon force = total reactive carbon force - coefficient * coeffcient * table value for reactive carbon and methanol - coefficient * coeffcient * table value for reactive carbon and methanol 
			if( r < 10.0){
			r0 = r/0.05;

			cgforce1 = linear_interp(r0,(r0+0.05),f[r0][1],f[r0+1][1],r);
			cgforce2 = linear_interp(r0,(r0+0.05),f[r0][1],f[r0+1][1],r);


			fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
			fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
			fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
			}
///////////////////////////////////////////

			//calculate distance between methanol and bonded chloride

			disx = x[i]-x[1003];
			disy = y[i]-y[1003];
			disz = z[i]-z[1003];

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

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

			//total bonded chlorine force = total bonded chlorine force - coefficient * coeffcient * table value for bonded chlorine and methanol - coeffcient *coeffcient * table value for reactive chlorine  and methanol
			if( r < 10.0){
			r0 = r/0.05;

			cgforce1 = linear_interp(r0,(r0+0.05),f[r0][2],f[r0+1][2],r);
			cgforce2 = linear_interp(r0,(r0+0.05),f[r0][0],f[r0+1][0],r);


			fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
			fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
			fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
			}
/////////////////////////////////////////////
			
			//methanol and propyl bead

			disx = x[i]-x[1004];
			disy = y[i]-y[1004];
			disz = z[i]-z[1004];

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

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

			//total propyl force = total propyl  force - coefficient * coeffcient * table value for propyl and methanol - coefficient * coeffcient * table value for propyl and methanol
			if( r < 10.0){
			r0 = r/0.05;

			cgforce1 = linear_interp(r0,(r0+0.05),f[r0][3],f[r0+1][3],r);
			cgforce2 = linear_interp(r0,(r0+0.05),f[r0][3],f[r0+1][3],r);


			fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
			fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
			fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
			}	
		}

///////////////////////////////////////////

		//positive ion and reactive chloride

		disx = x[1001]-x[1000];
		disy = y[1001]-y[1000];
		disz = z[1001]-z[1000];

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

//		fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

		if( r < 10.0){
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][4],f[r0+1][4],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][6],f[r0+1][6],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////

		//positive ion and reactive carbon bead

		disx = x[1001]-x[1002];
		disy = y[1001]-y[1002];
		disz = z[1001]-z[1002];

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

//		fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

		if( r < 10.0){
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][5],f[r0+1][5],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][5],f[r0+1][5],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////

		//positive ion and bonded chlorine

		disx = x[1001]-x[1003];
		disy = y[1001]-y[1003];
		disz = z[1001]-z[1003];

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

//		fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

		if( r < 10.0){
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][6],f[r0+1][6],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][4],f[r0+1][4],r);

		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////

		//positive ion and propyl bead

		disx = x[1001]-x[1004];
		disy = y[1001]-y[1004];
		disz = z[1001]-z[1004];

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

//		fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//		fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//		fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)


		if( r < 10.0){
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][7],f[r0+1][7],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][7],f[r0+1][7],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////
	
		//reactive chloride with reactive carbon bead

		disx = x[1000]-x[1002];
		disy = y[1000]-y[1002];
		disz = z[1000]-z[1002];
	
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
		r0 = r/0.05;
		r1 = r/0.01;
		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][8],f[r0+1][8],r);
		cgforce2 = linear_interp(r0,(r0+0.05),fbond[r1][1],fbond[r1+1][10],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);

		fx[1002] = fx[1002] - c[0]*c[0]*(cgforce1*disx/r)-c[1]*c[1]*(cgforce2*disx/r);
		fy[1002] = fy[1002] - c[0]*c[0]*(cgforce1*disy/r)-c[1]*c[1]*(cgforce2*disy/r);
		fz[1002] = fz[1002] - c[0]*c[0]*(cgforce1*disz/r)-c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////
	
		//reactive chloride with bonded chloride

		disx = x[1000]-x[1003];
		disy = y[1000]-y[1003];
		disz = z[1000]-z[1003];
	
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
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][9],f[r0+1][9],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][9],f[r0+1][9],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);

		fx[1002] = fx[1002] - c[0]*c[0]*(cgforce1*disx/r)-c[1]*c[1]*(cgforce2*disx/r);
		fy[1002] = fy[1002] - c[0]*c[0]*(cgforce1*disy/r)-c[1]*c[1]*(cgforce2*disy/r);
		fz[1002] = fz[1002] - c[0]*c[0]*(cgforce1*disz/r)-c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////
/*
		//reactive chloride with propyl bead

		disx = x[1000]-x[1004];
		disy = y[1000]-y[1004];
		disz = z[1000]-z[1004];
	
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
		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),f[r0][9],f[r0+1][9],r);



		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r);

		fx[1002] = fx[1002] - c[0]*c[0]*(cgforce1*disx/r);
		fy[1002] = fy[1002] - c[0]*c[0]*(cgforce1*disy/r);
		fz[1002] = fz[1002] - c[0]*c[0]*(cgforce1*disz/r);
		}
*/
///////////////////////////////////////////

		//reactive carbon and bonded chloride

		disx = x[1002]-x[1003];
		disy = y[1002]-y[1003];
		disz = z[1002]-z[1003];
	
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
		r0 = r/0.05;
		r1 = r/0.01;

		cgforce1 = linear_interp(r0,(r0+0.05),fbond[r1][1],fbond[r1+1][1],r);
		cgforce2 = linear_interp(r0,(r0+0.05),f[r0][8],f[r0+1][8],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);

		fx[1002] = fx[1002] - c[0]*c[0]*(cgforce1*disx/r)-c[1]*c[1]*(cgforce2*disx/r);
		fy[1002] = fy[1002] - c[0]*c[0]*(cgforce1*disy/r)-c[1]*c[1]*(cgforce2*disy/r);
		fz[1002] = fz[1002] - c[0]*c[0]*(cgforce1*disz/r)-c[1]*c[1]*(cgforce2*disz/r);
		}
///////////////////////////////////////////

		//reactive carbon and propyl

		disx = x[1002]-x[1004];
		disy = y[1002]-y[1004];
		disz = z[1002]-z[1004];
	
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

		r0 = r/0.05;

		cgforce1 = linear_interp(r0,(r0+0.05),fbond[r0][1],fbond[r0+1][1],r);
		cgforce2 = linear_interp(r0,(r0+0.05),fbond[r0][1],fbond[r0+1][1],r);


		fx[1000] = fx[1000] + c[0]*c[0]*(cgforce1*disx/r)+c[1]*c[1]*(cgforce2*disx/r);
		fy[1000] = fy[1000] + c[0]*c[0]*(cgforce1*disy/r)+c[1]*c[1]*(cgforce2*disy/r);
		fz[1000] = fz[1000] + c[0]*c[0]*(cgforce1*disz/r)+c[1]*c[1]*(cgforce2*disz/r);

		fx[1002] = fx[1002] - c[0]*c[0]*(cgforce1*disx/r)-c[1]*c[1]*(cgforce2*disx/r);
		fy[1002] = fy[1002] - c[0]*c[0]*(cgforce1*disy/r)-c[1]*c[1]*(cgforce2*disy/r);
		fz[1002] = fz[1002] - c[0]*c[0]*(cgforce1*disz/r)-c[1]*c[1]*(cgforce2*disz/r);

///////////////////////////////////////////
/*
		//bonded chloride and propyl

		disx = x[1003]-x[1004];
		disy = y[1003]-y[1004];
		disz = z[1003]-z[1004];
	
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

		r0=r/0.05;

		cgforce2 = linear_interp(r0,(r0+0.05),fbond[r0][1],fbond[r0+1][1],r)


		fx[1003] = fx[1003]-c[1]*c[1]*(cgforce2*disx/r);
		fy[1003] = fy[1003]-c[1]*c[1]*(cgforce2*disy/r);
		fz[1003] = fz[1003]-c[1]*c[1]*(cgforce2*disz/r);

		fx[1004] = fx[1004]+c[1]*c[1]*(cgforce2*disx/r);
		fy[1004] = fy[1004]+c[1]*c[1]*(cgforce2*disy/r);
		fz[1004] = fz[1004]+c[1]*c[1]*(cgforce2*disz/r);
		}
*/
///////////////////////////////////////////
/*

		//angle with reactive 

		disx = x[1000]-x[1002];
		disy = y[1000]-y[1002];
		disz = z[1000]-z[1002];

		disx2 = x[1004]-x[1002];
		disy2 = y[1004]-y[1002];
		disz2 = z[1004]-z[1002];		

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

		fx[1000] = fx[1000] - c[0]*c[0]*(a11*delx + a12*delx2);
		fy[1000] = fy[1000] - c[0]*c[0]*(a11*dely + a12*dely2);
		fz[1000] = fz[1000] - c[0]*c[0]*(a11*delz + a12*delz2);

		fx[1002] = fx[1002] + c[0]*c[0]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
		fy[1002] = fy[1002] + c[0]*c[0]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
		fz[1002] = fz[1002] + c[0]*c[0]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

		fx[1004] = fx[1004] - c[0]*c[0]*(a22*delx2 + a12*delx);
		fy[1004] = fy[1004] - c[0]*c[0]*(a22*dely2 + a12*dely);
		fz[1004] = fz[1004] - c[0]*c[0]*(a22*delz2 + a12*delz);

////////////////////////////////

		//angle with bonded

		disx = x[1003]-x[1002];
		disy = y[1003]-y[1002];
		disz = z[1003]-z[1002];

		disx2 = x[1004]-x[1002];
		disy2 = y[1004]-y[1002];
		disz2 = z[1004]-z[1002];		

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

		fx[1003] = fx[1003] - c[1]*c[1]*(a11*delx + a12*delx2);
		fy[1003] = fy[1003] - c[1]*c[1]*(a11*dely + a12*dely2);
		fz[1003] = fz[1003] - c[1]*c[1]*(a11*delz + a12*delz2);

		fx[1002] = fx[1002] + c[1]*c[1]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
		fy[1002] = fy[1002] + c[1]*c[1]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
		fz[1002] = fz[1002] + c[1]*c[1]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

		fx[1004] = fx[1004] - c[1]*c[1]*(a22*delx2 + a12*delx);
		fy[1004] = fy[1004] - c[1]*c[1]*(a22*dely2 + a12*dely);
		fz[1004] = fz[1004] - c[1]*c[1]*(a22*delz2 + a12*delz);
*/

		//weigh forces

		fx[1000] = fx[1000]*c[0]*c[1];
		fy[1000] = fy[1000]*c[0]*c[1];
		fz[1000] = fz[1000]*c[0]*c[1];

		fx[1002] = fx[1002]*c[0]*c[1];
		fy[1002] = fy[1002]*c[0]*c[1];
		fz[1002] = fz[1002]*c[0]*c[1];

		fx[1003] = fx[1003]*c[0]*c[1];
		fy[1003] = fy[1003]*c[0]*c[1];
		fz[1003] = fz[1003]*c[0]*c[1];
		
		fx[1004] = fx[1004]*c[0]*c[1];
		fy[1004] = fy[1004]*c[0]*c[1];
		fz[1004] = fz[1004]*c[0]*c[1];

	}

return 0;

}

double linear_interp(double x0,double x1,double y0,double y1,double x)
{
	//linearly interpolates between two points on a table
	double a = (y1-y0) / (x1-x0);
	double y = (y0 + (x-x0)*a);
	return y;
}

double nonreactive_radial_force(double c0,double c1,double r,double dis,double cgforce1,double cgforce2)
{

	//general weighter
	double newforce;
	newforce = c0*c0*(cgforce1*dis/r) + c1*c1*(cgforce2*dis/r);
	return newforce;
}
