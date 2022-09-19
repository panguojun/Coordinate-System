{#sim game gui and shader fake code
    params{param1:1;param2:2;param3:3}
    #gui
    gui{
        btn1{type:button;p:100,10;onclick:param1++}
        btn2{type:button;p:100,30;onclick:param2++}
    }
    #shader
    shader2d{
        shaderparams{mouse:1,1;param1:0;param2:0}
        polar_coord{params:float x, float y;
            return: coord(vec2(cos(x),sin(y)),vec2(sin(x),cos(y)));
        }
        mainimage{
            coord1:=polar_coord(x,y);   // coord
            dx:=x-0.25;dy:=y-0.5;   
            r1:=sqrt(dx*dx+dy*dy);
            dx=x-0.75;dy:y-0.5;
            r2:=sqrt(dx*dx+dy*dy);
            
            f:=sin(r1)+sin(r2);         // field
            color=rgb(f,f,0);           // display
        }
    }
}
