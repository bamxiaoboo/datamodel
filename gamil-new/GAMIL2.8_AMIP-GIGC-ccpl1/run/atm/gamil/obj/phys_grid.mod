	  x  æ   k820309    ë          13.0        `BZ                                                                                                           
       /data3/work/yuxinzhu/test/gamil-new/model_platform/models/atm/GAMIL2.0_128x60/src/physics/cam1/phys_grid.F90 PHYS_GRID              NPTHREADS CHUNKS NLCHUNKS LCHUNKS KNUHCS PGCOLS GS_COL_NUM GS_COL_OFFSET BTOFC_BLK_NUM BTOFC_CHK_NUM BTOFC_BLK_OFFSET BTOFC_CHK_OFFSET                                                     
       R8 SHR_KIND_R8 R4 SHR_KIND_R4                                                     
       PCOLS PVER BEGCHUNK ENDCHUNK                                                     
       PLON PLAT BEGLAT ENDLAT                      @                              
       NPES                      @                              
       LAT_LOCAL_PROC_ID                                                    
                        !                                                                                                      !                                                                                                      !                            	                                                      16                 !                            
                                                      26          @    !                                                  @    !                                                         !                                                                                  128                 !                                                                   <               60                 !                                                         !                                                    @   !                                                      !                                                               &                   &                                                                                                              #MPIPRIV1%MPI_BOTTOM    #MPIPRIV1%MPI_IN_PLACE    #MPIPRIV1%MPI_STATUS_IGNORE                                                                                                                                                                                                      p          p            p                                                                                                     #MPIPRIV2%MPI_STATUSES_IGNORE    #MPIPRIV2%MPI_ERRCODES_IGNORE                                                                              p          p          p            p          p                                                                                                           p          p            p                                                                                                     #MPIPRIVC%MPI_ARGVS_NULL    #MPIPRIVC%MPI_ARGV_NULL    -                                                                          p          p          p            p          p                                  -                                                                         p          p            p                                                    @                                'Ì                    #NCOLS    #LON    #GCOL     #LAT !   #OWNER "   #LCHUNK #                                                                                                                                             p          p            p                                                                                   D                   p          p            p                                                                      !                               p          p            p                                                                      "     Ä                                                         #     È                        @@                               $                       @@                               %                       @@                               &                              @                           '     '                    #CHUNKID (   #COL )                                               (                                                               )                                    @                           *     '                    #CHUNK +   #CCOL ,                                               +                                                               ,                             @@                               -                       @                                .                              @                           /     'h                    #NCOLS 0   #NLVLS 1   #PTER 2                                               0                                                               1                                                            2                                         &                   &                                                      @                                3                       @                                4                       @                               5     <              
      p          p <           p <                                    @                                6     <                    p          p <           p <                                    @                               7                   
      p          p          p <           p          p <                                   @                                8                      @                                9                     @                                :                                   &                                                    @                                ;                                   &                   &                                                    @                                <                                   &                   &                                                    @                                =                                   &                   &                                                    @                                >                                   &                                           #         @                                   ?                    #PHYS_GRID_INIT%ASSOCIATED @                                                                                                                                                                               @     ASSOCIATED #         @                                  A                   #CREATE_CHUNKS%MIN B   #CREATE_CHUNKS%MOD C   #OPT D   #CHUNKS_PER_THREAD E                                                                                                B     MIN                                            C     MOD           
                                  D                     
                                  E           #         @                                  F                   #ASSIGN_CHUNKS%MOD G   #OPT H                                                                      G     MOD           
                                  H           #         @                                   I                    #INDEX_BEG J   #INDEX_END K             D                                 J                      D                                 K            %         @                                L                           #LCHUNKID M             
                                  M           #         @                                   N                    #LCHUNKID O   #LATDIM P   #LATS Q             
                                  O                     
                                  P                    D                                 Q                         p          5  p        r P       5  p        r P                     #         @                                   R                    #LCHUNKID S   #LTH T   #COLS U   #LATS V             
                                  S                     
                                  T                    
                                  U                        p          5  p        r T       5  p        r T                              D                                 V                         p          5  p        r T       5  p        r T                     %         @                                W                           #LCHUNKID X   #COL Y             
                                  X                     
                                  Y           #         @                                   Z                    #LCHUNKID [   #LONDIM \   #LONS ]             
                                  [                     
                                  \                    D                                 ]                         p          5  p        r \       5  p        r \                     #         @                                   ^                    #LCHUNKID _   #LTH `   #COLS a   #LONS b             
                                  _                     
                                  `                    
                                  a                        p          5  p        r `       5  p        r `                              D                                 b                          p          5  p        r `       5  p        r `                     %         @                                c                           #LCHUNKID d   #COL e             
                                  d                     
                                  e           #         @                                   f                    #LCHUNKID g   #RLATDIM h   #RLATS i             
                                  g                     
                                  h                    D                                i                    
 !    p          5  p        r h       5  p        r h                     #         @                                   j                    #LCHUNKID k   #LTH l   #COLS m   #RLATS n             
                                  k                     
                                  l                    
                                  m                     "   p          5  p        r l       5  p        r l                              D                                n                    
 #    p          5  p        r l       5  p        r l                     %         @                                o                    
       #LCHUNKID p   #COL q             
                                  p                     
                                  q           #         @                                   r                    #LCHUNKID s   #RLONDIM t   #RLONS u             
                                  s                     
                                  t                    D                                u                    
 $    p          5  p        r t       5  p        r t                     #         @                                   v                    #LCHUNKID w   #LTH x   #COLS y   #RLONS z             
                                  w                     
                                  x                    
                                  y                     %   p          5  p        r x       5  p        r x                              D                                z                    
 &    p          5  p        r x       5  p        r x                     %         @                                {                    
       #LCHUNKID |   #COL }             
                                  |                     
                                  }           %         @                                ~                           #LONI    #LATJ              
                                                       
                                             %         @                                                           #IDX              
                                             #         @                                                       #LTH    #XYLONS    #XYLATS    #CKCOLS    #CKCIDS              
                                                      
                                                       '   p          5  p        r        5  p        r                               
                                                       (   p          5  p        r        5  p        r                               D                                                      )    p          5  p        r        5  p        r                               D                                                      *    p          5  p        r        5  p        r                      %         @                                                           #GCOL              
                                             #         @                                                       #LCID    #LATDIM    #GCOLS              
                                                       
                                                       D                                                     +              &                                           #         @                                                      #FDIM    #MDIM    #LDIM    #NLOND    #GLOBALFIELD    #LOCALCHUNKS              
                                                       
                                                       
                                                       
                                                      
                                                     
 ,           p        p <       p        5  p        r    p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r      p <         5  p        r        5  p        r      5  p        r      5  p        r      p <         5  p        r                               D                                                    
 -            p          5 r    5 r    p        5 r    5  p        r    p        p        p        5  p        r    p          5  p        r      p          5  p        r      & 5 r    5 r      5  p        r        5  p        r      p          5  p        r        5 r    5 r    p          5  p        r                      #         @                                                       #FDIM    #MDIM    #LDIM    #NLOND    #GLOBALFIELD    #LOCALCHUNKS              
                                                       
                                                       
                                                       
                                                      
                                                     	 2           p        p <       p        5  p        r    p        5  p        r    p        5  p        r    p          5  p        r      5  p        r      5  p        r      p <         5  p        r        5  p        r      5  p        r      5  p        r      p <         5  p        r                               D                                                    	 3            p          5 r    5 r    p        5 r    5  p        r    p        p        p        5  p        r    p          5  p        r      p          5  p        r      & 5 r    5 r      5  p        r        5  p        r      p          5  p        r        5 r    5 r    p          5  p        r                      #         @                                                       #FDIM    #MDIM    #LDIM     #NLOND ¡   #GLOBALFIELD ¢   #LOCALCHUNKS £             
                                                       
                                                       
                                                        
                                  ¡                    
                                  ¢                     8           p        p <       p        5  p        r    p        5  p        r ¡   p        5  p        r    p          5  p        r      5  p        r ¡     5  p        r      p <         5  p        r         5  p        r      5  p        r ¡     5  p        r      p <         5  p        r                                D                                 £                     9            p          5 r    5 r    p        5 r    5  p        r    p        p        p        5  p        r    p          5  p        r      p          5  p        r      & 5 r    5 r      5  p        r         5  p        r      p          5  p        r        5 r    5 r    p          5  p        r                       #         @                                  ¤                    #FDIM ¥   #MDIM ¦   #LDIM §   #NLOND ¨   #LOCALCHUNKS ©   #GLOBALFIELD ª             
                                  ¥                     
                                  ¦                     
                                  §                     
                                  ¨                    
                                 ©                    
 >           p          5 r    5 r    p        5 r    5  p        r ¦   p        p        p        5  p        r ¥   p          5  p        r ¥     p          5  p        r ¦     & 5 r    5 r      5  p        r §       5  p        r ¥     p          5  p        r ¦       5 r    5 r    p          5  p        r §                              D                                ª                    
 ?            p        p <       p        5  p        r ¦   p        5  p        r ¨   p        5  p        r ¥   p          5  p        r ¥     5  p        r ¨     5  p        r ¦     p <         5  p        r §       5  p        r ¥     5  p        r ¨     5  p        r ¦     p <         5  p        r §                     #         @                                   «                    #FDIM ¬   #MDIM ­   #LDIM ®   #NLOND ¯   #LOCALCHUNKS °   #GLOBALFIELD ±             
                                  ¬                     
                                  ­                     
                                  ®                     
                                  ¯                    
                                 °                    	 D           p          5 r    5 r    p        5 r    5  p        r ­   p        p        p        5  p        r ¬   p          5  p        r ¬     p          5  p        r ­     & 5 r    5 r      5  p        r ®       5  p        r ¬     p          5  p        r ­       5 r    5 r    p          5  p        r ®                              D                                ±                    	 E            p        p <       p        5  p        r ­   p        5  p        r ¯   p        5  p        r ¬   p          5  p        r ¬     5  p        r ¯     5  p        r ­     p <         5  p        r ®       5  p        r ¬     5  p        r ¯     5  p        r ­     p <         5  p        r ®                     #         @                                   ²                    #FDIM ³   #MDIM ´   #LDIM µ   #NLOND ¶   #LOCALCHUNKS ·   #GLOBALFIELD ¸             
                                  ³                     
                                  ´                     
                                  µ                     
                                  ¶                    
                                  ·                     J           p          5 r    5 r    p        5 r    5  p        r ´   p        p        p        5  p        r ³   p          5  p        r ³     p          5  p        r ´     & 5 r    5 r      5  p        r µ       5  p        r ³     p          5  p        r ´       5 r    5 r    p          5  p        r µ                              
D                                 ¸                     K            p        p <       p        5  p        r ´   p        5  p        r ¶   p        5  p        r ³   p          5  p        r ³     5  p        r ¶     5  p        r ´     p <         5  p        r µ       5  p        r ³     5  p        r ¶     5  p        r ´     p <         5  p        r µ                     #         @                                   ¹                    #IU º   #FDIM »   #MDIM ¼   #LDIM ½   #LOCALCHUNKS ¾             
                                  º                     
  @                               »                     
  @                               ¼                     
  @                               ½                    
  @                              ¾                    
 P           p          5 r    5 r    p        5 r    5  p        r ¼   p        p        p        5  p        r »   p          5  p        r »     p          5  p        r ¼     & 5 r    5 r      5  p        r ½       5  p        r »     p          5  p        r ¼       5 r    5 r    p          5  p        r ½                     #         @                                   ¿                    #IU À   #FDIM Á   #MDIM Â   #LDIM Ã   #LOCALCHUNKS Ä             
                                  À                     
  @                               Á                     
  @                               Â                     
  @                               Ã                    D @                              Ä                    
 R            p          5 r    5 r    p        5 r    5  p        r Â   p        p        p        5  p        r Á   p          5  p        r Á     p          5  p        r Â     & 5 r    5 r      5  p        r Ã       5  p        r Á     p          5  p        r Â       5 r    5 r    p          5  p        r Ã                     #         @                                   Å                    #RECORD_SIZE Æ   #BLOCK_BUFFER Ç   #CHUNK_BUFFER È             
                                  Æ                    
@ @                              Ç                    
 T   p           5  p        r Æ   5 r 3        5  p        r Æ   5 r 3                              D @                              È                    
 U    p           5  p        r Æ   5 r 4        5  p        r Æ   5 r 4                     #         @                                   É                    #BLOCKID Ê   #FDIM Ë   #LDIM Ì   #RECORD_SIZE Í   #PTER Î             
                                  Ê                     
                                  Ë                     
                                  Ì                     
                                  Í                    D                                 Î                     Z      p        5  p        r Ë   p          5  p        r Ë     5  p        r Ì       5  p        r Ë     5  p        r Ì                     #         @                                   Ï                    #LCHUNKID Ð   #FDIM Ñ   #LDIM Ò   #RECORD_SIZE Ó   #PTER Ô             
                                  Ð                     
                                  Ñ                     
                                  Ò                     
                                  Ó                    D                                 Ô                     [      p        5  p        r Ñ   p          5  p        r Ñ     5  p        r Ò       5  p        r Ñ     5  p        r Ò                     #         @                                   Õ                    #RECORD_SIZE Ö   #CHUNK_BUFFER ×   #BLOCK_BUFFER Ø             
                                  Ö                    D @                              ×                    
 \    p           5  p        r Ö   5 r 4        5  p        r Ö   5 r 4                              D @                              Ø                    
 ]    p           5  p        r Ö   5 r 3        5  p        r Ö   5 r 3                     #         @                                   Ù                    #LCHUNKID Ú   #FDIM Û   #LDIM Ü   #RECORD_SIZE Ý   #PTER Þ             
                                  Ú                     
                                  Û                     
                                  Ü                     
                                  Ý                    D                                 Þ                     b      p        5  p        r Û   p          5  p        r Û     5  p        r Ü       5  p        r Û     5  p        r Ü                     #         @                                   ß                    #BLOCKID à   #FDIM á   #LDIM â   #RECORD_SIZE ã   #PTER ä             
                                  à                     
                                  á                     
                                  â                     
                                  ã                    D                                 ä                     c      p        5  p        r á   p          5  p        r á     5  p        r â       5  p        r á     5  p        r â                                  fn#fn         b   uapp(PHYS_GRID    ¶  ^   J   SHR_KIND_MOD      ]   j  PPGRID    q  X   J  PMGRID    É  E   j  SPMD_DYN      R   J  MPI_GAMIL    `  @   J   MPISHORTHAND ,      p       R8+SHR_KIND_MOD=SHR_KIND_R8 ,     p       R4+SHR_KIND_MOD=SHR_KIND_R4      r       PCOLS+PPGRID    ò  r       PVER+PPGRID     d  @       BEGCHUNK+PPGRID     ¤  @       ENDCHUNK+PPGRID    ä  s       PLON+PMGRID    W  r       PLAT+PMGRID    É  @       BEGLAT+PMGRID    	  @       ENDLAT+PMGRID    I  @       NPES+SPMD_DYN ,     ¤       LAT_LOCAL_PROC_ID+MPI_GAMIL <   -  ¤      MPISHORTHAND!MPIPRIV1+MPISHORTHAND=MPIPRIV1 1   Ñ  H      MPIPRIV1%MPI_BOTTOM+MPISHORTHAND 3   	  H      MPIPRIV1%MPI_IN_PLACE+MPISHORTHAND 8   a	  ¤      MPIPRIV1%MPI_STATUS_IGNORE+MPISHORTHAND <   
        MPISHORTHAND!MPIPRIV2+MPISHORTHAND=MPIPRIV2 :   
  Ä      MPIPRIV2%MPI_STATUSES_IGNORE+MPISHORTHAND :   ]  ¤      MPIPRIV2%MPI_ERRCODES_IGNORE+MPISHORTHAND <           MPISHORTHAND!MPIPRIVC+MPISHORTHAND=MPIPRIVC 5     Ä      MPIPRIVC%MPI_ARGVS_NULL+MPISHORTHAND 4   N  ¤      MPIPRIVC%MPI_ARGV_NULL+MPISHORTHAND    ò         CHUNK      H   a   CHUNK%NCOLS    È     a   CHUNK%LON    d     a   CHUNK%GCOL          a   CHUNK%LAT      H   a   CHUNK%OWNER    ä  H   a   CHUNK%LCHUNK    ,  @       NLTHREADS    l  @       NGTHREADS    ¬  @       NCHUNKS    ì  f       KNUHC    R  H   a   KNUHC%CHUNKID      H   a   KNUHC%COL    â  e       COLUMN_MAP !   G  H   a   COLUMN_MAP%CHUNK       H   a   COLUMN_MAP%CCOL    ×  @       NGCOLS      @       NLCOLS    W  p       BTOFC_PTERS "   Ç  H   a   BTOFC_PTERS%NCOLS "     H   a   BTOFC_PTERS%NLVLS !   W  ¬   a   BTOFC_PTERS%PTER       @       BLOCK_BUF_NRECS     C  @       CHUNK_BUF_NRECS             CLAT_P             NLON_P    «  ´       CLON_P    _  @       PHYSGRID_SET      @       LOCAL_DP_MAP    ß         CHUNK_NCOLS    k  ¤       CHUNK_LAT      ¤       CHUNK_LON    ³  ¤       CHUNK_PID     W         LCHUNK_TO_CHUNK    ã  è       PHYS_GRID_INIT *   Ë  C      PHYS_GRID_INIT%ASSOCIATED      È       CREATE_CHUNKS "   Ö  <      CREATE_CHUNKS%MIN "     <      CREATE_CHUNKS%MOD "   N  @   a   CREATE_CHUNKS%OPT 0     @   a   CREATE_CHUNKS%CHUNKS_PER_THREAD    Î         ASSIGN_CHUNKS "   N  <      ASSIGN_CHUNKS%MOD "     @   a   ASSIGN_CHUNKS%OPT $   Ê  f       GET_CHUNK_INDICES_P .   0   @   a   GET_CHUNK_INDICES_P%INDEX_BEG .   p   @   a   GET_CHUNK_INDICES_P%INDEX_END    °   ^       GET_NCOLS_P %   !  @   a   GET_NCOLS_P%LCHUNKID    N!  l       GET_LAT_ALL_P '   º!  @   a   GET_LAT_ALL_P%LCHUNKID %   ú!  @   a   GET_LAT_ALL_P%LATDIM #   :"  ´   a   GET_LAT_ALL_P%LATS    î"  s       GET_LAT_VEC_P '   a#  @   a   GET_LAT_VEC_P%LCHUNKID "   ¡#  @   a   GET_LAT_VEC_P%LTH #   á#  ´   a   GET_LAT_VEC_P%COLS #   $  ´   a   GET_LAT_VEC_P%LATS    I%  g       GET_LAT_P #   °%  @   a   GET_LAT_P%LCHUNKID    ð%  @   a   GET_LAT_P%COL    0&  l       GET_LON_ALL_P '   &  @   a   GET_LON_ALL_P%LCHUNKID %   Ü&  @   a   GET_LON_ALL_P%LONDIM #   '  ´   a   GET_LON_ALL_P%LONS    Ð'  s       GET_LON_VEC_P '   C(  @   a   GET_LON_VEC_P%LCHUNKID "   (  @   a   GET_LON_VEC_P%LTH #   Ã(  ´   a   GET_LON_VEC_P%COLS #   w)  ´   a   GET_LON_VEC_P%LONS    +*  g       GET_LON_P #   *  @   a   GET_LON_P%LCHUNKID    Ò*  @   a   GET_LON_P%COL    +  n       GET_RLAT_ALL_P (   +  @   a   GET_RLAT_ALL_P%LCHUNKID '   À+  @   a   GET_RLAT_ALL_P%RLATDIM %    ,  ´   a   GET_RLAT_ALL_P%RLATS    ´,  t       GET_RLAT_VEC_P (   (-  @   a   GET_RLAT_VEC_P%LCHUNKID #   h-  @   a   GET_RLAT_VEC_P%LTH $   ¨-  ´   a   GET_RLAT_VEC_P%COLS %   \.  ´   a   GET_RLAT_VEC_P%RLATS    /  g       GET_RLAT_P $   w/  @   a   GET_RLAT_P%LCHUNKID    ·/  @   a   GET_RLAT_P%COL    ÷/  n       GET_RLON_ALL_P (   e0  @   a   GET_RLON_ALL_P%LCHUNKID '   ¥0  @   a   GET_RLON_ALL_P%RLONDIM %   å0  ´   a   GET_RLON_ALL_P%RLONS    1  t       GET_RLON_VEC_P (   2  @   a   GET_RLON_VEC_P%LCHUNKID #   M2  @   a   GET_RLON_VEC_P%LTH $   2  ´   a   GET_RLON_VEC_P%COLS %   A3  ´   a   GET_RLON_VEC_P%RLONS    õ3  g       GET_RLON_P $   \4  @   a   GET_RLON_P%LCHUNKID    4  @   a   GET_RLON_P%COL "   Ü4  d       GET_CHUNK_OWNER_P '   @5  @   a   GET_CHUNK_OWNER_P%LONI '   5  @   a   GET_CHUNK_OWNER_P%LATJ    À5  Y       CHUNK_INDEX     6  @   a   CHUNK_INDEX%IDX "   Y6         GET_CHUNK_COORD_P &   Ú6  @   a   GET_CHUNK_COORD_P%LTH )   7  ´   a   GET_CHUNK_COORD_P%XYLONS )   Î7  ´   a   GET_CHUNK_COORD_P%XYLATS )   8  ´   a   GET_CHUNK_COORD_P%CKCOLS )   69  ´   a   GET_CHUNK_COORD_P%CKCIDS !   ê9  Z       GET_GCOL_OWNER_P &   D:  @   a   GET_GCOL_OWNER_P%GCOL    :  i       GET_GCOL_ALL_P $   í:  @   a   GET_GCOL_ALL_P%LCID &   -;  @   a   GET_GCOL_ALL_P%LATDIM %   m;     a   GET_GCOL_ALL_P%GCOLS '   ù;         SCATTER_FIELD_TO_CHUNK ,   <  @   a   SCATTER_FIELD_TO_CHUNK%FDIM ,   Ì<  @   a   SCATTER_FIELD_TO_CHUNK%MDIM ,   =  @   a   SCATTER_FIELD_TO_CHUNK%LDIM -   L=  @   a   SCATTER_FIELD_TO_CHUNK%NLOND 3   =  D  a   SCATTER_FIELD_TO_CHUNK%GLOBALFIELD 3   Ð?  d  a   SCATTER_FIELD_TO_CHUNK%LOCALCHUNKS (   4B         SCATTER_FIELD_TO_CHUNK4 -   ÇB  @   a   SCATTER_FIELD_TO_CHUNK4%FDIM -   C  @   a   SCATTER_FIELD_TO_CHUNK4%MDIM -   GC  @   a   SCATTER_FIELD_TO_CHUNK4%LDIM .   C  @   a   SCATTER_FIELD_TO_CHUNK4%NLOND 4   ÇC  D  a   SCATTER_FIELD_TO_CHUNK4%GLOBALFIELD 4   F  d  a   SCATTER_FIELD_TO_CHUNK4%LOCALCHUNKS +   oH         SCATTER_FIELD_TO_CHUNK_INT 0   I  @   a   SCATTER_FIELD_TO_CHUNK_INT%FDIM 0   BI  @   a   SCATTER_FIELD_TO_CHUNK_INT%MDIM 0   I  @   a   SCATTER_FIELD_TO_CHUNK_INT%LDIM 1   ÂI  @   a   SCATTER_FIELD_TO_CHUNK_INT%NLOND 7   J  D  a   SCATTER_FIELD_TO_CHUNK_INT%GLOBALFIELD 7   FL  d  a   SCATTER_FIELD_TO_CHUNK_INT%LOCALCHUNKS &   ªN         GATHER_CHUNK_TO_FIELD +   =O  @   a   GATHER_CHUNK_TO_FIELD%FDIM +   }O  @   a   GATHER_CHUNK_TO_FIELD%MDIM +   ½O  @   a   GATHER_CHUNK_TO_FIELD%LDIM ,   ýO  @   a   GATHER_CHUNK_TO_FIELD%NLOND 2   =P  d  a   GATHER_CHUNK_TO_FIELD%LOCALCHUNKS 2   ¡R  D  a   GATHER_CHUNK_TO_FIELD%GLOBALFIELD '   åT         GATHER_CHUNK_TO_FIELD4 ,   xU  @   a   GATHER_CHUNK_TO_FIELD4%FDIM ,   ¸U  @   a   GATHER_CHUNK_TO_FIELD4%MDIM ,   øU  @   a   GATHER_CHUNK_TO_FIELD4%LDIM -   8V  @   a   GATHER_CHUNK_TO_FIELD4%NLOND 3   xV  d  a   GATHER_CHUNK_TO_FIELD4%LOCALCHUNKS 3   ÜX  D  a   GATHER_CHUNK_TO_FIELD4%GLOBALFIELD *    [         GATHER_CHUNK_TO_FIELD_INT /   ³[  @   a   GATHER_CHUNK_TO_FIELD_INT%FDIM /   ó[  @   a   GATHER_CHUNK_TO_FIELD_INT%MDIM /   3\  @   a   GATHER_CHUNK_TO_FIELD_INT%LDIM 0   s\  @   a   GATHER_CHUNK_TO_FIELD_INT%NLOND 6   ³\  d  a   GATHER_CHUNK_TO_FIELD_INT%LOCALCHUNKS 6   _  D  a   GATHER_CHUNK_TO_FIELD_INT%GLOBALFIELD '   [a         WRITE_FIELD_FROM_CHUNK *   Úa  @   a   WRITE_FIELD_FROM_CHUNK%IU ,   b  @   a   WRITE_FIELD_FROM_CHUNK%FDIM ,   Zb  @   a   WRITE_FIELD_FROM_CHUNK%MDIM ,   b  @   a   WRITE_FIELD_FROM_CHUNK%LDIM 3   Úb  d  a   WRITE_FIELD_FROM_CHUNK%LOCALCHUNKS &   >e         READ_CHUNK_FROM_FIELD )   ½e  @   a   READ_CHUNK_FROM_FIELD%IU +   ýe  @   a   READ_CHUNK_FROM_FIELD%FDIM +   =f  @   a   READ_CHUNK_FROM_FIELD%MDIM +   }f  @   a   READ_CHUNK_FROM_FIELD%LDIM 2   ½f  d  a   READ_CHUNK_FROM_FIELD%LOCALCHUNKS )   !i  }       TRANSPOSE_BLOCK_TO_CHUNK 5   i  @   a   TRANSPOSE_BLOCK_TO_CHUNK%RECORD_SIZE 6   Þi  Ô   a   TRANSPOSE_BLOCK_TO_CHUNK%BLOCK_BUFFER 6   ²j  Ô   a   TRANSPOSE_BLOCK_TO_CHUNK%CHUNK_BUFFER *   k         BLOCK_TO_CHUNK_SEND_PTERS 2   
l  @   a   BLOCK_TO_CHUNK_SEND_PTERS%BLOCKID /   Jl  @   a   BLOCK_TO_CHUNK_SEND_PTERS%FDIM /   l  @   a   BLOCK_TO_CHUNK_SEND_PTERS%LDIM 6   Êl  @   a   BLOCK_TO_CHUNK_SEND_PTERS%RECORD_SIZE /   
m  $  a   BLOCK_TO_CHUNK_SEND_PTERS%PTER *   .n         BLOCK_TO_CHUNK_RECV_PTERS 3   ³n  @   a   BLOCK_TO_CHUNK_RECV_PTERS%LCHUNKID /   ón  @   a   BLOCK_TO_CHUNK_RECV_PTERS%FDIM /   3o  @   a   BLOCK_TO_CHUNK_RECV_PTERS%LDIM 6   so  @   a   BLOCK_TO_CHUNK_RECV_PTERS%RECORD_SIZE /   ³o  $  a   BLOCK_TO_CHUNK_RECV_PTERS%PTER )   ×p  }       TRANSPOSE_CHUNK_TO_BLOCK 5   Tq  @   a   TRANSPOSE_CHUNK_TO_BLOCK%RECORD_SIZE 6   q  Ô   a   TRANSPOSE_CHUNK_TO_BLOCK%CHUNK_BUFFER 6   hr  Ô   a   TRANSPOSE_CHUNK_TO_BLOCK%BLOCK_BUFFER *   <s         CHUNK_TO_BLOCK_SEND_PTERS 3   Ás  @   a   CHUNK_TO_BLOCK_SEND_PTERS%LCHUNKID /   t  @   a   CHUNK_TO_BLOCK_SEND_PTERS%FDIM /   At  @   a   CHUNK_TO_BLOCK_SEND_PTERS%LDIM 6   t  @   a   CHUNK_TO_BLOCK_SEND_PTERS%RECORD_SIZE /   Át  $  a   CHUNK_TO_BLOCK_SEND_PTERS%PTER *   åu         CHUNK_TO_BLOCK_RECV_PTERS 2   iv  @   a   CHUNK_TO_BLOCK_RECV_PTERS%BLOCKID /   ©v  @   a   CHUNK_TO_BLOCK_RECV_PTERS%FDIM /   év  @   a   CHUNK_TO_BLOCK_RECV_PTERS%LDIM 6   )w  @   a   CHUNK_TO_BLOCK_RECV_PTERS%RECORD_SIZE /   iw  $  a   CHUNK_TO_BLOCK_RECV_PTERS%PTER 