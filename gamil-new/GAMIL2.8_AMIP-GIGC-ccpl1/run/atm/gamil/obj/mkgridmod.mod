	  �)  J   k820309    �          13.0        ^BZ                                                                                                           
       /data3/work/yuxinzhu/test/gamil-new/model_platform/models/lnd/CLM2/src/mksrfdata/mkgridMod.F90 MKGRIDMOD          @       � @                              
       R8 SHR_KIND_R8                                                    
                                                          
                                                          
                            @                              
       GETFIL                                                    
                                                           
       MASTERPROC                  !                                                                                     #         @       !                            	                   #GETFIL%LEN_TRIM 
   #GETFIL%TRIM    #GETFIL%PRESENT    #FULPATH    #LOCFN    #IFLAG                  @                            
     LEN_TRIM               @                                 TRIM               @                                 PRESENT           
   @                                                 1             @                                                  1           
  @                                                          !                                                                                                                            <               60           @                                                        @ @                                    <                    p          p <           p <                                   @                                                                                                                                �               128           @ @                                                 
      p �         p �         p <           p �         p <                                    @ @                                                 
      p �         p �         p <           p �         p <                                    @                                                          p �         p �         p <           p �         p <                                    @                                                   
      p �         p �         p <           p �         p <                         #         @                                                       #NLAT    #NLON    #NUMLON    #LONGXY    #LATIXY    #EDGEN     #EDGEE !   #EDGES "   #EDGEW #   #LATS $   #LONW %             
                                                       
                                                      
                                                       7   p          5 O p            5 O p                                   
                                                     
 8     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                   
                                                     
 9     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                    
                                       
                
                                 !     
                
                                 "     
                
                                 #     
                                               $                    
 :    p           5 O p        n                                       1     5 O p        n                                      1                                                                    %                    
 ;      p         5 O p        n                                           1p           5 O p        n                                      1  5 O p             5 O p        n                                      1  5 O p                                    #         @                                  &                    #NLAT '   #NLON (   #NUMLON )   #LONGXY *   #LATIXY +   #LATS ,   #LONW -             
                                  '                     
                                  (                    
                                  )                     <   p          5 O p            5 O p                                   
                                 *                    
 =     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                   
                                 +                    
 >     p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                                                   ,                    
 ?    p           5 O p        n                                       1     5 O p        n                                      1                                                                    -                    
 @      p         5 O p        n                                           1p           5 O p        n                                      1  5 O p             5 O p        n                                      1  5 O p                                    #         @                                   .                
   #CELLAREA_REGIONAL%ABS /   #CELLAREA_REGIONAL%SIN 0   #NLAT 1   #NLON 2   #NUMLON 3   #LATS 4   #LONW 5   #EDGEN 6   #EDGEE 7   #EDGES 8   #EDGEW 9   #AREA :                                               /     ABS                                             0     SIN           
                                  1                     
                                  2                    
                                  3                     /   p          5 O p            5 O p                                   
                                 4                    
 0   p           5 O p        n                                       1     5 O p        n                                      1                                    
                                 5                    
 1     p         5 O p        n                                           1p           5 O p        n                                      1  5 O p             5 O p        n                                      1  5 O p                                              
                                 6     
                
                                 7     
                
                                 8     
                
                                 9     
                                               :                    
 2      p        5 O p        p          5 O p          5 O p            5 O p          5 O p                          #         @                                  ;                   #CELLAREA_GLOBAL%SIN <   #NLAT =   #NLON >   #NUMLON ?   #LATS @   #LONW A   #AREA B                                               <     SIN           
                                  =                     
                                  >                    
                                  ?                     3   p          5 O p            5 O p                                   
                                 @                    
 4   p           5 O p        n                                       1     5 O p        n                                      1                                    
                                 A                    
 5     p         5 O p        n                                           1p           5 O p        n                                      1  5 O p             5 O p        n                                      1  5 O p                                                                             B                    
 6      p        5 O p        p          5 O p          5 O p            5 O p          5 O p                          #         @                                   C                   #MKGRID_CAM%ABS D   #CAM_LONGXY E   #CAM_LATIXY F   #CAM_NUMLON G   #CAM_LANDFRAC H   #CAM_LANDMASK I                                              D     ABS           
                                 E                   
              &                   &                                                     
                                 F                   
              &                   &                                                     
                                  G                                 &                                                     
                                 H                   
              &                   &                                                     
                                  I                                 &                   &                                              �   q      fn#fn      O   J   SHR_KIND_MOD    `  @   J   CLM_VARPAR    �  @   J   CLM_VARSUR    �  @   J   CLM_VARCTL       G   J  FILEUTILS    g  @   J   AREAMOD    �  K   J  SPMDMOD ,   �  p       R8+SHR_KIND_MOD=SHR_KIND_R8 !   b  �       GETFIL+FILEUTILS 3     A      GETFIL%LEN_TRIM+FILEUTILS=LEN_TRIM +   H  =      GETFIL%TRIM+FILEUTILS=TRIM 1   �  @      GETFIL%PRESENT+FILEUTILS=PRESENT )   �  L   e   GETFIL%FULPATH+FILEUTILS '     L   e   GETFIL%LOCFN+FILEUTILS '   ]  @   e   GETFIL%IFLAG+FILEUTILS "   �  @       MASTERPROC+PMGRID "   �  r       LSMLAT+CLM_VARPAR '   O  @       POLE_POINTS+CLM_VARSUR "   �  �       NUMLON+CLM_VARSUR $   #  @       FULLGRID+CLM_VARSUR "   c  s       LSMLON+CLM_VARPAR "   �  �       LONGXY+CLM_VARSUR "   �  �       LATIXY+CLM_VARSUR $   >	  �       LANDMASK+CLM_VARSUR $   �	  �       LANDFRAC+CLM_VARSUR *   �
  �       CELLEDGE_REGIONAL+AREAMOD /   f  @   a   CELLEDGE_REGIONAL%NLAT+AREAMOD /   �  @   a   CELLEDGE_REGIONAL%NLON+AREAMOD 1   �  �   a   CELLEDGE_REGIONAL%NUMLON+AREAMOD 1   �  �   a   CELLEDGE_REGIONAL%LONGXY+AREAMOD 1   �  �   a   CELLEDGE_REGIONAL%LATIXY+AREAMOD 0   �  @   a   CELLEDGE_REGIONAL%EDGEN+AREAMOD 0   �  @   a   CELLEDGE_REGIONAL%EDGEE+AREAMOD 0     @   a   CELLEDGE_REGIONAL%EDGES+AREAMOD 0   B  @   a   CELLEDGE_REGIONAL%EDGEW+AREAMOD /   �    a   CELLEDGE_REGIONAL%LATS+AREAMOD /   �  �  a   CELLEDGE_REGIONAL%LONW+AREAMOD (   ?  �       CELLEDGE_GLOBAL+AREAMOD -   �  @   a   CELLEDGE_GLOBAL%NLAT+AREAMOD -     @   a   CELLEDGE_GLOBAL%NLON+AREAMOD /   S  �   a   CELLEDGE_GLOBAL%NUMLON+AREAMOD /   �  �   a   CELLEDGE_GLOBAL%LONGXY+AREAMOD /   �  �   a   CELLEDGE_GLOBAL%LATIXY+AREAMOD -   �    a   CELLEDGE_GLOBAL%LATS+AREAMOD -     �  a   CELLEDGE_GLOBAL%LONW+AREAMOD *   �  �       CELLAREA_REGIONAL+AREAMOD .   �  <      CELLAREA_REGIONAL%ABS+AREAMOD .   �  <      CELLAREA_REGIONAL%SIN+AREAMOD /     @   a   CELLAREA_REGIONAL%NLAT+AREAMOD /   L  @   a   CELLAREA_REGIONAL%NLON+AREAMOD 1   �  �   a   CELLAREA_REGIONAL%NUMLON+AREAMOD /   0    a   CELLAREA_REGIONAL%LATS+AREAMOD /   F  �  a   CELLAREA_REGIONAL%LONW+AREAMOD 0   �  @   a   CELLAREA_REGIONAL%EDGEN+AREAMOD 0   -  @   a   CELLAREA_REGIONAL%EDGEE+AREAMOD 0   m  @   a   CELLAREA_REGIONAL%EDGES+AREAMOD 0   �  @   a   CELLAREA_REGIONAL%EDGEW+AREAMOD /   �  �   a   CELLAREA_REGIONAL%AREA+AREAMOD (   �  �       CELLAREA_GLOBAL+AREAMOD ,   �   <      CELLAREA_GLOBAL%SIN+AREAMOD -   �   @   a   CELLAREA_GLOBAL%NLAT+AREAMOD -   !  @   a   CELLAREA_GLOBAL%NLON+AREAMOD /   D!  �   a   CELLAREA_GLOBAL%NUMLON+AREAMOD -   �!    a   CELLAREA_GLOBAL%LATS+AREAMOD -   �"  �  a   CELLAREA_GLOBAL%LONW+AREAMOD -   �$  �   a   CELLAREA_GLOBAL%AREA+AREAMOD    �%  �       MKGRID_CAM    Q&  <      MKGRID_CAM%ABS &   �&  �   a   MKGRID_CAM%CAM_LONGXY &   1'  �   a   MKGRID_CAM%CAM_LATIXY &   �'  �   a   MKGRID_CAM%CAM_NUMLON (   a(  �   a   MKGRID_CAM%CAM_LANDFRAC (   )  �   a   MKGRID_CAM%CAM_LANDMASK 