function F = F(oxygen_file,hydrogen_file)

O = load(oxygen_file);
H = load(hydrogen_file);
T = length(O);

RO = O(:,1:3);
RH = H(:,1:3);
ROH = RO - RH;

mO = 16;
mH = 1;
mu = mO*mH/(mO + mH);

FO = O(:,4:6);
FH = H(:,4:6);
dF = (FO/mO - FH/mH);

for t = 1:T
    x = 0;
    for k = 1:3
        x = x + ROH(t,k)*FH(t,k);
    end
    F(t) = x*mu;
end

end