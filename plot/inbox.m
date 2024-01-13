function [stns,stns_id]=inbox(stns1,stn_rg)
j=1;
for i=1:length(stns1)
    stns0=stns1(i,:);
    if stns0(1)>stn_rg(1) && stns0(1)<stn_rg(2) && stns0(2)>stn_rg(3) && stns0(2)<stn_rg(4)
        stns(j,:)=stns0;
        stns_id(j)=i;
        j=j+1;
    end
end