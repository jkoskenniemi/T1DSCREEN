# T1DSCREEN

A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr

OPEN IN EDIT MODE TO DISPLAY THE FILE WORKFLOW CORRECTLY

Original files───────────────────>Gather.R──────────────>Harmonize.R───────────────────────>/analyses/*.Rmd		

assocvariants.annotated.txt		
	      ├─────────────────────────IL6ST_anno──┐		
	      ├─────────────────────────IL2RA_anno──│┐		
	      ├─────────────────────────IL2RB_anno──││┐		
	      ├─────────────────────────IL2RG_anno──│││┐		
      	      ├─────────────────────────CXCL10_anno─││││┐		
	      ├─────────────────────────TYK2_anno───│││││┐		
	      ├─────────────────────────TNF_anno────││││││┐		
      	      ├─────────────────────────JAK2_anno───│││││││┐	
	      ├─────────────────────────IL12B_anno──││││││││┐	
	      └─────────────────────────IL6R_anno───│││││││││┐	
					            ││││││││││
Ferkingstad_data			            ││││││││││	
  	2620_4_IL6ST_gp130__soluble.txt─IL6ST_prot──┴│││││││││─IL6ST_prot_anno──┐
 	 3151_6_IL2RA_IL_2_sRa.txt───────IL2RA_prot───┴││││││││─IL2RA_prot_anno──│┐
 	 5260_80_TYK2_TYK2.txt───────────IL2RB_prot────┴│││││││─IL2RB_prot_anno──││┐	
 	 2634_2_IL2RG_IL_2_sRg.txt───────IL2RG_prot─────┴││││││─IL2RG_prot_anno──│││┐	
	 4141_79_CXCL10_IP_10.txt────────CXCL10_prot─────┴│││││─CXCL10_prot_anno─││││┐		
	 5260_80_TYK2_TYK2.txt───────────TYK2_prot────────┴││││─TYK2_prot_anno───│││││┐		
	 5936_53_TNF_TNF_a.txt───────────TNF_prot──────────┴│││─TNF_prot_anno────││││││┐	
 	 11816_84_JAK2_JAK2.txt──────────JAK2_prot──────────┴││─JAK_prot_anno────│││││││┐ 	
	 13733_5_IL12B_IL_12_p40.txt─────IL12B_prot──────────┴│─IL12B_prot_anno──││││││││┐ 
       			  		  		                         │││││││││ 
Chiou et al. 2022 Nature				       		         │││││││││ 
34012112-GCST90014023-EFO_0001359.h.tsv"				         │││││││││ 
				       ├──IL6ST_T1D──────────────────────────────┴││││││││──IL6ST───>IL6ST.rmd	
                        	       ├──IL2RA_T1D───────────────────────────────┴│││││││──IL2RA───>IL2RA.rmd	
				       ├──IL2RB_T1D────────────────────────────────┴││││││──IL2RB───>IL2RB.rmd
	                               ├──IL2RG_T1D─────────────────────────────────┴│││││──IL2RG───>IL2RG.rmd
				       ├──CXCL10_T1D─────────────────────────────────┴││││──CXCL10──>CXCL10.rmd
                        	       ├──TYK2_T1D────────────────────────────────────┴│││──TYK2────>TYK2.rmd  
                        	       ├──TNF_T1D──────────────────────────────────────┴││──TNF─────>TNF.rmd
                        	       ├──JAK2_T1D──────────────────────────────────────┴│──JAK2────>JAK2.rmd
                        	       ├──IL12B_T1D──────────────────────────────────────┴──IL12B───>IL12B.rmd
                        	       └──IL6R_T1D───────────────────────────────────────┬──IL6R────>CRP.rmd
	                                                 				 │
35459240-GCST90029070-EFO_0004458.h.tsv"                                 		 │
         15602_43_IL6R_IL_6_sRa.txt──────IL6R_crp──────────────IL6R_crp_anno─────────────┘
							      
