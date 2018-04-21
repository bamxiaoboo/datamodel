! Declarations

      integer ntim
      parameter(ntim=453)
      integer yrdata(ntim)      ! yearly data values
      real(r8)    sconst(ntim)      ! input time-varying solar const (W/m2)

! Ouput type declaration

      character*64, parameter :: ramp_type = &
                       'RAMP_SCON: using ramp_scon data'
      logical :: ramp_write
      data ramp_write / .true. /

! Input data values

      data yrdata / 1849  ,&
           1850  ,1851  ,1852  ,1853  ,1854  , &
           1855  ,1856  ,1857  ,1858  ,1859  , &
           1860  ,1861  ,1862  ,1863  ,1864  , &
           1865  ,1866  ,1867  ,1868  ,1869  , & 
           1870  ,1871  ,1872  ,1873  ,1874  , &
           1875  ,1876  ,1877  ,1878  ,1879  , &
           1880  ,1881  ,1882  ,1883  ,1884  , &
           1885  ,1886  ,1887  ,1888  ,1889  , &
           1890  ,1891  ,1892  ,1893  ,1894  , &
           1895  ,1896  ,1897  ,1898  ,1899  , &
           1900  ,1901  ,1902  ,1903  ,1904  , &
           1905  ,1906  ,1907  ,1908  ,1909  , &
           1910  ,1911  ,1912  ,1913  ,1914  , &
           1915  ,1916  ,1917  ,1918  ,1919  , &
           1920  ,1921  ,1922  ,1923  ,1924  , &
           1925  ,1926  ,1927  ,1928  ,1929  , &
           1930  ,1931  ,1932  ,1933  ,1934  , &
           1935  ,1936  ,1937  ,1938  ,1939  , &
           1940  ,1941  ,1942  ,1943  ,1944  , &
           1945  ,1946  ,1947  ,1948  ,1949  , &
           1950  ,1951  ,1952  ,1953  ,1954  , &
           1955  ,1956  ,1957  ,1958  ,1959  , &
           1960  ,1961  ,1962  ,1963  ,1964  , &
           1965  ,1966  ,1967  ,1968  ,1969  , &
           1970  ,1971  ,1972  ,1973  ,1974  , &
           1975  ,1976  ,1977  ,1978  ,1979  , &
           1980  ,1981  ,1982  ,1983  ,1984  , &
           1985  ,1986  ,1987  ,1988  ,1989  , &
           1990  ,1991  ,1992  ,1993  ,1994  , &
           1995  ,1996  ,1997  ,1998  ,1999  , &
           2000  ,2001  ,2002  ,2003  ,2004  , &
           2005  ,2006  ,2007  ,2008  ,2009  , &
           2010  ,2011  ,2012  ,2013  ,2014  , &
           2015  ,2016  ,2017  ,2018  ,2019  , &
           2020  ,2021  ,2022  ,2023  ,2024  , &
           2025  ,2026  ,2027  ,2028  ,2029  , &
           2030  ,2031  ,2032  ,2033  ,2034  , &
           2035  ,2036  ,2037  ,2038  ,2039  , &
           2040  ,2041  ,2042  ,2043  ,2044  , &
           2045  ,2046  ,2047  ,2048  ,2049  , &
           2050  ,2051  ,2052  ,2053  ,2054  , &
           2055  ,2056  ,2057  ,2058  ,2059  , &
           2060  ,2061  ,2062  ,2063  ,2064  , &
           2065  ,2066  ,2067  ,2068  ,2069  , &
           2070  ,2071  ,2072  ,2073  ,2074  , &
           2075  ,2076  ,2077  ,2078  ,2079  , &
           2080  ,2081  ,2082  ,2083  ,2084  , &
           2085  ,2086  ,2087  ,2088  ,2089  , &
           2090  ,2091  ,2092  ,2093  ,2094  , &
           2095  ,2096  ,2097  ,2098  ,2099  , &
           2100  ,2101  ,2102  ,2103  ,2104  , &
           2105  ,2106  ,2107  ,2108  ,2109  , &
           2110  ,2111  ,2112  ,2113  ,2114  , &
           2115  ,2116  ,2117  ,2118  ,2119  , &
           2120  ,2121  ,2122  ,2123  ,2124  , &
           2125  ,2126  ,2127  ,2128  ,2129  , &
           2130  ,2131  ,2132  ,2133  ,2134  , &
           2135  ,2136  ,2137  ,2138  ,2139  , &
           2140  ,2141  ,2142  ,2143  ,2144  , &
           2145  ,2146  ,2147  ,2148  ,2149  , &
           2150  ,2151  ,2152  ,2153  ,2154  , &
           2155  ,2156  ,2157  ,2158  ,2159  , &
           2160  ,2161  ,2162  ,2163  ,2164  , &
           2165  ,2166  ,2167  ,2168  ,2169  , &
           2170  ,2171  ,2172  ,2173  ,2174  , &
           2175  ,2176  ,2177  ,2178  ,2179  , &
           2180  ,2181  ,2182  ,2183  ,2184  , &
           2185  ,2186  ,2187  ,2188  ,2189  , &
           2190  ,2191  ,2192  ,2193  ,2194  , &
           2195  ,2196  ,2197  ,2198  ,2199  , &
           2200  ,2201  ,2202  ,2203  ,2204  , &
           2205  ,2206  ,2207  ,2208  ,2209  , &
           2210  ,2211  ,2212  ,2213  ,2214  , &
           2215  ,2216  ,2217  ,2218  ,2219  , &
           2220  ,2221  ,2222  ,2223  ,2224  , &
           2225  ,2226  ,2227  ,2228  ,2229  , &
           2230  ,2231  ,2232  ,2233  ,2234  , &
           2235  ,2236  ,2237  ,2238  ,2239  , &
           2240  ,2241  ,2242  ,2243  ,2244  , &
           2245  ,2246  ,2247  ,2248  ,2249  , &
           2250  ,2251  ,2252  ,2253  ,2254  , &
           2255  ,2256  ,2257  ,2258  ,2259  , &
           2260  ,2261  ,2262  ,2263  ,2264  , &
           2265  ,2266  ,2267  ,2268  ,2269  , &
           2270  ,2271  ,2272  ,2273  ,2274  , &
           2275  ,2276  ,2277  ,2278  ,2279  , &
           2280  ,2281  ,2282  ,2283  ,2284  , &
           2285  ,2286  ,2287  ,2288  ,2289  , &
           2290  ,2291  ,2292  ,2293  ,2294  , &
           2295  ,2296  ,2297  ,2298  ,2299  , &
           2300  ,2301  / 

! solar constanst 1849-2009. Use 1997 value for 1997-2100

       data  sconst /    1.36594409e+06,  &
        1.36568229e+06,  1.36571609e+06,  1.36565933e+06,  1.36557026e+06,  1.36546148e+06,   &
        1.36533702e+06,  1.36442273e+06,  1.36399050e+06,  1.36496666e+06,  1.36562630e+06,   &
        1.36585932e+06,  1.36577695e+06,  1.36537661e+06,  1.36537835e+06,  1.36545350e+06,   &
        1.36542752e+06,  1.36538891e+06,  1.36535975e+06,  1.36551936e+06,  1.36573903e+06,   &
        1.36596940e+06,  1.36591968e+06,  1.36585429e+06,  1.36561579e+06,  1.36551237e+06,   &
        1.36537619e+06,  1.36520717e+06,  1.36523048e+06,  1.36525175e+06,  1.36530540e+06,   &
        1.36547495e+06,  1.36562976e+06,  1.36566759e+06,  1.36453338e+06,  1.36204571e+06,   &
        1.36391950e+06,  1.36443543e+06,  1.36436693e+06,  1.36471700e+06,  1.36452393e+06,   &
        1.36429409e+06,  1.36476617e+06,  1.36508359e+06,  1.36552673e+06,  1.36577209e+06,   &
        1.36570771e+06,  1.36510192e+06,  1.36499395e+06,  1.36509200e+06,  1.36522768e+06,   &
        1.36526782e+06,  1.36521046e+06,  1.36473572e+06,  1.36366252e+06,  1.36492641e+06,   &
        1.36521657e+06,  1.36554845e+06,  1.36538265e+06,  1.36542273e+06,  1.36546624e+06,   &
        1.36535383e+06,  1.36530979e+06,  1.36481775e+06,  1.36471959e+06,  1.36514056e+06,   &
        1.36559315e+06,  1.36583447e+06,  1.36599379e+06,  1.36596391e+06,  1.36573222e+06,   &
        1.36538661e+06,  1.36532520e+06,  1.36532392e+06,  1.36539948e+06,  1.36537486e+06,   &
        1.36559346e+06,  1.36570499e+06,  1.36591130e+06,  1.36569534e+06,  1.36553750e+06,   &
        1.36561340e+06,  1.36552896e+06,  1.36533012e+06,  1.36525307e+06,  1.36543376e+06,   &
        1.36564186e+06,  1.36605621e+06,  1.36600094e+06,  1.36588225e+06,  1.36588685e+06,   &
        1.36584498e+06,  1.36579721e+06,  1.36563785e+06,  1.36547575e+06,  1.36555978e+06,   &
        1.36582892e+06,  1.36593957e+06,  1.36616435e+06,  1.36630943e+06,  1.36617836e+06,   &
        1.36594434e+06,  1.36573061e+06,  1.36567596e+06,  1.36554919e+06,  1.36556675e+06,   &
        1.36573564e+06,  1.36628380e+06,  1.36666389e+06,  1.36663286e+06,  1.36638285e+06,   &
        1.36616629e+06,  1.36565114e+06,  1.36542134e+06,  1.36455252e+06,  1.36390599e+06,   &
        1.36465913e+06,  1.36533873e+06,  1.36575183e+06,  1.36548113e+06,  1.36538917e+06,   &
        1.36583009e+06,  1.36579761e+06,  1.36595666e+06,  1.36560742e+06,  1.36541250e+06,   &
        1.36478456e+06,  1.36530995e+06,  1.36570814e+06,  1.36608928e+06,  1.36640569e+06,   &
        1.36653365e+06,  1.36657424e+06,  1.36497346e+06,  1.36432006e+06,  1.36506034e+06,   &
        1.36532905e+06,  1.36530038e+06,  1.36553363e+06,  1.36589510e+06,  1.36649787e+06,   &
        1.36640328e+06,  1.36509981e+06,  1.36327499e+06,  1.36480975e+06,  1.36530335e+06,   &
        1.36545791e+06,  1.36545364e+06,  1.36561276e+06,  1.36603531e+06,  1.36636007e+06,   &
        1.36664035e+06,  1.36657720e+06,  1.36665570e+06,  1.36620500e+06,  1.36602300e+06,   &
        1.36582950e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561198e+06,   &
        1.36573985e+06,  1.36610198e+06,  1.36638507e+06,  1.36666535e+06,  1.36660220e+06,   &
        1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,   &
        1.36572400e+06,  1.36569180e+06,  1.36561198e+06,  1.36573985e+06,  1.36610198e+06,   &
        1.36638507e+06,  1.36666535e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,   &
        1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,   &
        1.36561198e+06,  1.36573985e+06,  1.36610198e+06,  1.36638507e+06,  1.36666535e+06,   &
        1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,   &
        1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561198e+06,  1.36573985e+06,   &
        1.36610198e+06,  1.36638507e+06,  1.36666535e+06,  1.36660220e+06,  1.36668070e+06,   &
        1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,   &
        1.36569180e+06,  1.36561198e+06,  1.36573985e+06,  1.36610198e+06,  1.36638507e+06,   &
        1.36666535e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,   &
        1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561198e+06,   &
        1.36573985e+06,  1.36610198e+06,  1.36638507e+06,  1.36666535e+06,  1.36660220e+06,   &
        1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,   &
        1.36572400e+06,  1.36569180e+06,  1.36561198e+06,  1.36573985e+06,  1.36610198e+06,   &
        1.36638507e+06,  1.36666535e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,   &
        1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,   &        
        1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,   &
        1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,   &
        1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,   &
        1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,   &
        1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,   &
        1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,   &
        1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,   &
        1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,   &
        1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,   &
        1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,   &
        1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,   &
        1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,   &
        1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,   &
        1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,   &
        1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,   &
        1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,   &
        1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,   &
        1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,   &
        1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,   &
        1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,   &
        1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,   &
        1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,   &
        1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,   &
        1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,   &
        1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,   &
        1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,   &
        1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,   &
        1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,   &
        1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,   &
        1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,   &
        1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,   &
        1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,   &
        1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,  1.36604800e+06,   &
        1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,  1.36561210e+06,   &
        1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,  1.36660220e+06,   &
        1.36668070e+06,  1.36623000e+06,  1.36604800e+06,  1.36585450e+06,  1.36581070e+06,   &
        1.36572400e+06,  1.36569180e+06,  1.36561210e+06,  1.36573990e+06,  1.36610210e+06,   &
        1.36638510e+06,  1.36668360e+06,  1.36660220e+06,  1.36668070e+06,  1.36623000e+06,   &
        1.36604800e+06,  1.36585450e+06,  1.36581070e+06,  1.36572400e+06,  1.36569180e+06,   &
        1.36561210e+06,  1.36573990e+06,  1.36610210e+06,  1.36638510e+06,  1.36668360e+06,   &
        1.36660220e+06,  1.36668070e+06 /