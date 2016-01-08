#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

			//total reactive chloride force = total reactive choride force - coefficient * coeffcient * table value for reactive chloride and methanol - coeffcient *coeffcient * table value for bonded chloride and methanol

			fx[1000] = fx[1000] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][]*disx/r)
			fy[1000] = fy[1000] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][]*disy/r)
			fz[1000] = fz[1000] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][]*disz/r)

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

//			fx[i] = fx[i] - c[t][0]*(f[r][0]*disx/r)-c[t][1](f[r][0]*disx/r)
//			fy[i] = fy[i] - c[t][0]*(f[r][0]*disy/r)-c[t][1](f[r][0]*disy/r)
//			fz[i] = fz[i] - c[t][0]*(f[r][0]*disz/r)-c[t][1](f[r][0]*disz/r)

			//total reactive carbon force = total reactive carbon force - coefficient * coeffcient * table value for reactive carbon and methanol - coefficient * coeffcient * table value for reactive carbon and methanol 

			fx[1002] = fx[1002] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
			fy[1002] = fy[1002] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
			fz[1002] = fz[1002] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

			fx[1003] = fx[1003] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
			fy[1003] = fy[1003] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
			fz[1003] = fz[1003] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

			fx[1004] = fx[1004] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
			fy[1004] = fy[1004] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
			fz[1004] = fz[1004] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)	
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



		fx[1000] = fx[1000] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1000] = fy[1000] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1000] = fz[1000] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

		fx[1002] = fx[1002] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1002] = fy[1002] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1002] = fz[1002] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

		fx[1003] = fx[1003] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[t][1]*(f[r][0]*disx/r)
		fy[1003] = fy[1003] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[t][1]*(f[r][0]*disy/r)
		fz[1003] = fz[1003] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[t][1]*(f[r][0]*disz/r)

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

		fx[1004] = fx[1004] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1004] = fy[1004] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1004] = fz[1004] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

		fx[1000] = fx[1000] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1000] = fy[1000] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1000] = fz[1000] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1002] = fx[1002] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[t][1]*(f[r][0]*disx/r)
		fy[1002] = fy[1002] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[t][1]*(f[r][0]*disy/r)
		fz[1002] = fz[1002] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[t][1]*(f[r][0]*disz/r)

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

		fx[1000] = fx[1000] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1000] = fy[1000] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1000] = fz[1000] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1003] = fx[1003] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1003] = fy[1003] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1003] = fz[1003] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

///////////////////////////////////////////

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

		fx[1000] = fx[1000] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1000] = fy[1000] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1000] = fz[1000] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1004] = fx[1004] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[t][1]*(f[r][0]*disx/r)
		fy[1004] = fy[1004] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[t][1]*(f[r][0]*disy/r)
		fz[1004] = fz[1004] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[t][1]*(f[r][0]*disz/r)

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

		fx[1002] = fx[1002] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1002] = fy[1002] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1002] = fz[1002] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1003] = fx[1003] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1003] = fy[1003] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1003] = fz[1003] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

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

		fx[1002] = fx[1002] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1002] = fy[1002] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1002] = fz[1002] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1004] = fx[1004] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1004] = fy[1004] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1004] = fz[1004] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

///////////////////////////////////////////

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

		fx[1003] = fx[1003] - c[0]*c[0]*(f[r][0]*disx/r)-c[1]*c[1]*(f[r][0]*disx/r)
		fy[1003] = fy[1003] - c[0]*c[0]*(f[r][0]*disy/r)-c[1]*c[1]*(f[r][0]*disy/r)
		fz[1003] = fz[1003] - c[0]*c[0]*(f[r][0]*disz/r)-c[1]*c[1]*(f[r][0]*disz/r)

		fx[1004] = fx[1004] + c[0]*c[0]*(f[r][0]*disx/r)+c[1]*c[1]*(f[r][0]*disx/r)
		fy[1004] = fy[1004] + c[0]*c[0]*(f[r][0]*disy/r)+c[1]*c[1]*(f[r][0]*disy/r)
		fz[1004] = fz[1004] + c[0]*c[0]*(f[r][0]*disz/r)+c[1]*c[1]*(f[r][0]*disz/r)

///////////////////////////////////////////
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

		fx[1000] = fx[1000] - c[0]*(a11*delx + a12*delx2);
		fy[1000] = fy[1000] - c[0]*(a11*dely + a12*dely2);
		fz[1000] = fz[1000] - c[0]*(a11*delz + a12*delz2);

		fx[1002] = fx[1002] + c[0]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
		fy[1002] = fy[1002] + c[0]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
		fz[1002] = fz[1002] + c[0]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

		fx[1004] = fx[1004] - c[0]*(a22*delx2 + a12*delx);
		fy[1004] = fy[1004] - c[0]*(a22*dely2 + a12*dely);
		fz[1004] = fz[1004] - c[0]*(a22*delz2 + a12*delz);

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

		fx[1003] = fx[1003] - c[1]*(a11*delx + a12*delx2);
		fy[1003] = fy[1003] - c[1]*(a11*dely + a12*dely2);
		fz[1003] = fz[1003] - c[1]*(a11*delz + a12*delz2);

		fx[1002] = fx[1002] + c[1]*(a11*delx + a12*delx2 + a22*delx2 + a12*delx);
		fy[1002] = fy[1002] + c[1]*(a11*dely + a12*dely2 + a22*dely2 + a12*dely);
		fz[1002] = fz[1002] + c[1]*(a11*delz + a12*delz2 + a22*delz2 + a12*delz);

		fx[1004] = fx[1004] - c[1]*(a22*delx2 + a12*delx);
		fy[1004] = fy[1004] - c[1]*(a22*dely2 + a12*dely);
		fz[1004] = fz[1004] - c[1]*(a22*delz2 + a12*delz);


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

