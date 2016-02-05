                     /************************************************/
                     /* PUSH !!! Check on angle to use neighbourhood */
                     /************************************************/
                     /****************************/
                     /* Grain Orientation Factor */
                     /* Neigbourhood Flipping    */
                     /****************************/
		              if (oriented_on)
		              {
       		           prob = drand48();
		                 if (prob < (bp->gr[*ngr]->ph))
		                 {
			                 /*fprintf(stderr,"H");*/
			                 tmp_n_neigh = 8;
#ifdef NASTAC_RULE
			                 tmp_oni = (prob < 0.5) ? bp->nbhd.onhl : bp->nbhd.onhr ;
#else /*RCA_RULE*/
			                 tmp_oni = (bp->gr[*ngr]->dir_angle > 0) ? bp->nbhd.onhl : bp->nbhd.onhr ;
#endif /*NASTAC_RULE*/
		                 }
		                 else if (prob < (bp->gr[*ngr]->ph+bp->gr[*ngr]->pq))
		                 {
			                 /*fprintf(stderr,"Q");*/
			                 tmp_n_neigh = 6;
			                 tmp_oni = bp->nbhd.onq;
		                 }
		                 else
		                 {
		                    /*fprintf(stderr,"8");*/
			                 tmp_n_neigh = 10;
			                 tmp_oni = bp->nbhd.ono;
		                 }
                      /* end neighbourhood flipping */
		              } /*end if oriented on*/
