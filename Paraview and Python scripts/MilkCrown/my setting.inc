#include "textures.inc"     
#include "colors.inc"
#include "functions.inc"     
#include "table.inc"   
#include "metals.inc"
#include "lamp.inc"

#declare milk= material {
texture {
    pigment { color rgbt <0.8627, 0.851, 0.8039, 1> }
    finish {
      ambient 0.0
      diffuse 0.0
      reflection {
        0.0, 1.0
        //metallic
      }

      specular 0.4
      roughness 0.002
    }
  }
  interior {
    ior 1.5
    dispersion 1.01
    fade_distance 0.1
    fade_power 2
    media {
        scattering { 
            2 <0.5, 0.65, 0.4> 
            extinction  1.0
        }       
    }   
  }     
}
#declare cam=camera {	
	location <0,-25,0.5>
	angle 35.0
	look_at <0,0,-1.5>  
	rotate<0,0,30>
}    
#declare ls=light_source{<70,70,70> color White}      
#declare ss=sky_sphere{
 pigment{ gradient <0,1,0>
          color_map{   
            [0.05 color White  ]
            [0.1 color Gray50 ] 
            [0.15 color Gray75 ]
            [0.2 color Black  ]
                   } // end color_map
     scale 3
     rotate <6,30,0>
     translate <60,60,0>
 } 
} 
#declare fg=fog{fog_type   2   distance 65  color rgb<1,0.99,0.9>
    fog_offset 0.1 fog_alt  2.0 turbulence 0.2}    
       

//===================================================================
    
    
global_settings {
	ambient_light color rgb <1.0, 1.0, 1.0>
	assumed_gamma 1.5    
	photons { count 5000 }
}   

light_source { 
	< 2.5*20, 100, 2.5*10 >
	color rgb < 1, 1, 1 > 
}

// sky sphere

sphere {
	<0,0,0>, 10000

	hollow

	pigment {
		color rgb 1
	}
}

camera {cam}    
//light_source{ls}        
//sky_sphere{ss} 
//fog{fg}           

#declare tank= 
difference{              
box{<-3.05,-3.05,-6>,<3.05,3.05,-3>}
box{<-3,-3,-6>,<3,3,4.1>}  
texture{Glass3}
}       

//=================================   


object{tank} 
 
object{
 TABLE
 translate <0,-2.4,0> 
 rotate<90,0,0>    
 scale 10
}

object{
FLOOR
rotate<90,0,0>
translate<0,0,-10>
}         

object{
CHAIR
rotate<90,0,-30>
translate<-3,2.8,-3>  
scale 5
} 

object{
Lamp_On          
rotate<90,0,0>
translate <0,2.5,-1>
scale 10
}
