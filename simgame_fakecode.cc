{   #gui
    gui{
        params{param1:1;param2:2;param3:3}
        btn1{type:button;p:100,10;onclick:onbtn1}
        btn2{type:button;p:100,30;onclick:btn2}
    }
    #shader
    shader2d{
        params{mouse:1,1}
        mainimage{
            dx1:x-0.25;dy:y-0.5;
            r1:sqrt(dx1*dx1+dy1*dy1);
            dx2:x-0.75;dy:y-0.5;
            r2:sqrt(dx2*dx2+dy2*dy2);
            f:sin(r1)+sin(r2);
            color: rgb(f,f,0);
        }
    }
}
