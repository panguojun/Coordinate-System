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
            coord1:polar_coord(x,y);
            dx1:x-0.25;dy:y-0.5;
            r1:sqrt(dx1*dx1+dy1*dy1);
            dx2:x-0.75;dy:y-0.5;
            r2:sqrt(dx2*dx2+dy2*dy2);
            f:sin(r1)+sin(r2);
            color: rgb(f,f,0);
        }
    }
}
