/************************************************/
/* PULL METHOD                                  */
/************************************************/
/* Run through all cells updating as needed.    */
/************************************************/
   for(k=0;k<nz;k++) {          /* loop cells in z direction */
      for(j=0;j<ny;j++) {       /* loop cells in y direction */
         for(i=0;i<nx;i++) {    /* loop cells in x direction */

            if ( *np == SOLID) {
	       /* do nothing for now... */

            } else if ( *np  == LIQUID) {
	       /* Check the neighbours, if solid, turn it! */
	       for (n=0; n<n_neigh; n++) {
	          if ( *(op+oni[n]) >= SOLID ) {
	            /* if ( *(op+oni[n]) >= (bp->gr[*(ogr + oni[n])].mdir[n] * SOLID ) ) { */
		     *np = EMBRYO;
		     *ngr = *(ogr + oni[n]);
		    /* } */
		  }
		}
              /* Continuous casting mode: calc. nucleation for every cell */
	      if (cell_nuc_on) {
		 Tunder = Tliq -  cell_temp_calc_cc (bp, sbnum, i, j);

                 if ( Tunder > 0) {
                   if (cell_nuc(bp, Tunder) == 1) { /* New Cell nucleated!!! */
		     if (bp->nprops.ngr == bp->nprops.gd_max_total) {
			fprintf(stderr,"ERROR: sb_ca_step - Max # grains set by user [%d] exceeded. Increase option MaxTotGrains.\n",bp->nprops.gd_max_total);
		     } else {
		        *np = EMBRYO;
                        bp->nprops.ngr +=1;
                        bp->sb[sbnum]->ngr +=1;
		        *ngr = bp->nprops.ngr;
                        /* init the grain structure for that grain... */
                        init_new_grain (bp, bp->nprops.ngr, sbnum, i+j*nx+k*nx*ny,1 );

	             }
	           }
	         }
	       }


            } else { /* Growing!!! */

              /* Continuous casting mode: calc. temp. for every cell */
	      if (cell_temp_on) {
		 Tunder = Tliq -  cell_temp_calc_cc (bp, sbnum, i, j);
                 if ( Tunder > 0) {
                    dr = growth(bp->mprops.gg_const,bp->mprops.gg_cub, bp->delt, Tunder);
                    dfs = dr / bp->size_c[0];
	         } else {
                    dr = 0.0;
                    dfs = 0.0;
	         }
	      }

              /* Grain Orientation Factor */
              dfs_tmp = oriented_on ?
                    (dfs * (bp->gr[*ngr].dir_fact)) : dfs;
	       /* add in the extra fs... */
	       *np  += dfs_tmp;

	       if ( *np  > 1.0 ) {
	           sum_del_fs -= *np - 1.0; /* only add fraction of dfs */
	           *np  = 1.0;
		   sp->ncsolid++;
	          /* TO DO: put in push... */
	       }
	       sum_del_fs += dfs;
            } /* end of check cell state */
	    np++;
	    op++;
	    ngr++; ogr++;
         }
	 op += 2;
	 ogr += 2;
      }
      op += 2*(nx+2);
      ogr += 2*(nx+2);

   } /* end of loop through all cells for pull method*/
