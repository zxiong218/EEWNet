function [delta, az, baz]=distaz(stalon,stalat,evtlon,evtlat)
        lat1=stalat;
        lon1=stalon;
        lat2=evtlat;
        lon2=evtlon;
        if (lat1 == lat2) && (lon1 == lon2)
            delta = 0.0;
            az = 0.0;
            baz = 0.0;
            return
        end
        
        rad=2.*pi/360.0;
%{
	c
	c scolat and ecolat are the geocentric colatitudes
	c as defined by Richter (pg. 318)
	c
	c Earth Flattening of 1/298.257 take from Bott (pg. 3)
	c
 %}
        sph=1.0/298.257;

        scolat=pi/2.0 - atan((1.-sph)*(1.-sph)*tan(lat1*rad));
        ecolat=pi/2.0 - atan((1.-sph)*(1.-sph)*tan(lat2*rad));
        slon=lon1*rad;
        elon=lon2*rad;
%{
	c
	c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
	c     These are defined for the pt. 1
	c
%}
        a=sin(scolat)*cos(slon);
        b=sin(scolat)*sin(slon);
        c=cos(scolat);
        d=sin(slon);
        e=-cos(slon);
        g=-c*e;
        h=c*d;
        k=-sin(scolat);
%{
	c
	c  aa - ee are the same as a - e, except for pt. 2
	c
%}
        aa=sin(ecolat)*cos(elon);
        bb=sin(ecolat)*sin(elon);
        cc=cos(ecolat);
        dd=sin(elon);
        ee=-cos(elon);
        gg=-cc*ee;
        hh=cc*dd;
        kk=-sin(ecolat);
%{
	c
	c  Bullen, Sec 10.2, eqn. 4
	c
%}
        delrad=acos(a*aa + b*bb + c*cc);
        delta=delrad/rad;
%{
	c
	c  Bullen, Sec 10.2, eqn 7 / eqn 8
	c
	c    pt. 1 is unprimed, so this is technically the baz
	c
	c  Calculate baz this way to avoid quadrant problems
	c
%}
        rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.;
        rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.;
        dbaz=atan2(rhs1,rhs2);
        if (dbaz<0.0)
            dbaz=dbaz+2*pi;
        end
        
        baz=dbaz/rad;
%{
	c
	c  Bullen, Sec 10.2, eqn 7 / eqn 8
	c
	c    pt. 2 is unprimed, so this is technically the az
	c
%}
        rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.;
        rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.;
        daz=atan2(rhs1,rhs2);
        if daz<0.0
	        daz=daz+2*pi;
        end
        
        az=daz/rad;
%{
	c
	c   Make sure 0.0 is always 0.0, not 360.
	c
%}
        if (abs(baz-360.) < .00001)
	        baz=0.0;
        end
        if (abs(az-360.) < .00001)
	        az=0.0;
        end
end