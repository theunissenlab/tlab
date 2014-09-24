function all_birds=getAllBird;

all_birds=[];

for sex=['m','f'];
    gender={};
    gender.sex=sex;
    gender.birds=[]
    for birdn=1:16;
        bird=getOneBird(sex,birdn);
        gender.birds=[gender.birds;bird];
    end;
    all_birds=[all_birds;gender];
end;
    