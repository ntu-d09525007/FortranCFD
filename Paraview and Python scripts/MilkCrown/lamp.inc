
//==========================================================             
// Lamp
//==========================================================             
// Textures
//==========================================================

#declare Light=light_source {<0,0,0>
    color rgb<1,0.8,0.6>*15
        area_light 0.1*x,y*0.1 8, 8 jitter adaptive 1
    fade_distance 0.3
    fade_power 2
}
             
#declare T_Shade_On=texture{
    pigment{
        bozo turbulence 0.6
        color_map{
           [0 rgbf<1,0.8,0,0.3>]
           [1 rgbf<1,0.8,0,0.7>]
        }
    }
    finish{ambient 0 diffuse 0.8 specular 0.05 roughness 0.01}
}
#declare T_Shade_Off=texture{
    pigment{
        bozo turbulence 0.6
        color_map{
           [0 rgb<1,0.6,0>]
           [1 rgb<1,0.8,0>]
        }
    }
    finish{ambient 0 diffuse 0.8 specular 0.05 roughness 0.01}
}
#declare T_Iron=texture{
    pigment{rgb <0.8,0.8,0.9>*0.5} 
    normal{agate bump_size 0.3 scale 0.01}
    finish{metallic brilliance 3 ambient 0 diffuse 1 specular 1 roughness 0.05
        reflection{0,0.2}
    }
}
#declare T_Douille=texture{T_Brass_3A finish{ambient 0}}

//==========================================================             
// CSG
//==========================================================             
#declare Bulb=sphere{0,1
    hollow
    texture{pigment{Clear}finish{ambient 0 diffuse 0}}
    interior{
        media{     
            emission White*10
            density{
                spherical
                poly_wave 5
                color_map{
                    [0 color Black]
                    [0.1 color White]
                }       
            }
            
        }
    }                         
    scale 0.1
}
#declare ElFootLamp=prism{
	cubic_spline
	0,
	1.000,
	47,
	<0.199,0.097>,	<0.186,0.087>,	<0.163,0.086>,	<0.144,0.086>,	<0.124,0.086>,
	<0.104,0.089>,	<0.084,0.086>,	<0.051,0.079>,	<0.019,0.065>,	<0.006,0.041>,
	<0.017,0.017>,	<0.047,0.004>,	<0.050,0.017>,	<0.039,0.020>,	<0.024,0.027>,
	<0.016,0.040>,	<0.025,0.058>,	<0.042,0.065>,	<0.063,0.071>,	<0.084,0.074>,
	<0.102,0.076>,	<0.125,0.075>,	<0.146,0.074>,	<0.180,0.072>,	<0.200,0.079>,
	<0.210,0.090>,	<0.217,0.109>,	<0.211,0.131>,	<0.202,0.140>,	<0.190,0.147>,
	<0.170,0.153>,	<0.148,0.152>,	<0.132,0.146>,	<0.123,0.138>,	<0.119,0.120>,
	<0.129,0.105>,	<0.149,0.098>,	<0.133,0.114>,	<0.132,0.131>,	<0.144,0.139>,
	<0.160,0.142>,	<0.182,0.137>,	<0.199,0.124>,	<0.203,0.111>,	<0.199,0.097>,
	<0.186,0.087>,	<0.163,0.086>
	scale <1,0.02,-1>
	translate z*0.153
	rotate -x*90
	translate z*0.01
	scale 0.8                                             
	translate x*0.01
}
#declare ColumnLamp=isosurface{
        #declare rPost=0.03;
        #declare yPost=0.7;
        function{f_helix1(x,y,z,1,100,0.35*rPost,0.55*rPost,2,1,0)}
        contained_by {box {<-rPost,0,-rPost>,<rPost,yPost,rPost>}}
}                  
#local DecoLeaves=intersection{
        difference{sphere{0,1} sphere{0,0.99}}
        plane{z,0 rotate y*20 inverse}        
        plane{z,0 rotate -y*20}        
        translate <-1,0.7,0> rotate z*20
}                   
#local RingLeaves=union{  
        #local teta=0;                                        
        #while (teta<360)
                object{DecoLeaves rotate z*50 translate -x*0.8 rotate y*teta}
                #local teta=teta+36;
        #end
        scale <0.06,0.06,0.06>*0.8
}                
                            
#declare JointLampBottom=union{
        object{RingLeaves scale <0.8,1,0.8> translate y*0.05}
        torus{0.04,0.01}
        sphere{0,0.041}
        cylinder{0,7*y,3.1 scale 1/100}
}
#declare FootLamp=union{
        box{<-0.07,-0.01,-0.07>,<0.07,0,0.07> translate y*0.05 rotate y*45}
        object{ElFootLamp}
        object{ElFootLamp rotate y*90}
        object{ElFootLamp rotate y*180}
        object{ElFootLamp rotate y*270}
        scale <1,1.3,1>
}

#declare JointLampTop=union{
        torus{0.034,0.008 translate y*0.085}
        #local teta=0;                                        
        #while (teta<360)
                torus{0.02,0.003 rotate x*90 translate <-0.02,0,0> rotate y*teta translate y*0.085}
                #local teta=teta+360/15;
        #end
        lathe{
        	cubic_spline
        	16,
        	<0.022,0.000>,	<0.067,0.000>,	<0.060,0.020>,	<0.067,0.057>,
        	<0.058,0.065>,	<0.066,0.074>,	<0.050,0.078>,	<0.053,0.101>,
        	<0.064,0.104>,	<0.054,0.119>,	<0.056,0.139>,	<0.045,0.148>,
        	<0.054,0.167>,	<0.047,0.185>,	<0.027,0.198>,	<0.009,0.201>
        	translate -y*0.198
        	rotate x*180
        	scale 0.5
        }
}                             


#declare rShadeDown=0.27*0.7;
#declare rShadeUp=0.15*0.7;
#declare yShade=0.23;
#declare rShadeRim=0.005;
#declare eShade=rShadeRim;
#declare Shade_On=union{
    difference{
        cone{0,rShadeDown,y*yShade,rShadeUp}
        cone{-0.0001*y,rShadeDown-eShade,y*(yShade+0.001),rShadeUp-eShade}
        texture{T_Shade_On scale 0.05}
        interior{
            media{
                scattering {1,2}
            }
        }
        hollow
    }                       
    union{
        torus{rShadeDown,rShadeRim}
        torus{rShadeUp,rShadeRim translate y*yShade}
        #declare i=0;
        #while (i<4)
            cylinder{-z*rShadeDown,y*yShade/10,rShadeRim rotate y*360*i/4}
            #declare i=i+1;
        #end
        texture{T_Brass_2A}
//                no_shadow
    }
//    rotate -z*3
}                       
        
#declare Shade_Off=union{
    difference{
        cone{0,rShadeDown,y*yShade,rShadeUp}
        cone{-0.0001*y,rShadeDown-eShade,y*(yShade+0.001),rShadeUp-eShade}
        texture{T_Shade_Off scale 0.05}
    }                       
    union{
        torus{rShadeDown,rShadeRim}
        torus{rShadeUp,rShadeRim translate y*yShade}
        #declare i=0;
        #while (i<4)
            cylinder{-z*rShadeDown,y*yShade/10,rShadeRim rotate y*360*i/4}
            #declare i=i+1;
        #end
        texture{T_Brass_2A}
//                no_shadow
    }
//    rotate -z*3
}                       
#declare Lamp_On=union{                     
    object{FootLamp} // h=0.07
    object{ColumnLamp scale <0.8,1,0.8> translate y*0.07} // h=0.7
    union{
        object{JointLampBottom translate y*0.77} // h=0.076
        object{JointLampTop texture{T_Douille} scale 1.2 translate y*0.816} // h=0.099
        scale <0.5,1,0.5>
//          no_shadow
    }
    union{
        object{Light translate y*0.16}
        object{Bulb translate y*0.08}
        object{Shade_On}
        translate y*0.9
    }
    texture{T_Iron}
    scale 1.5
}