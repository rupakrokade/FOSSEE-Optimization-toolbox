
c = [-750 -1000];
G = [
	1 1;
	1 2;
	4 3;
    ];
h = [10 15 25]';
A = [
	0.5 1
];
b=[7.5];
l = "a";
q = [];
e = []
dims=list("l",l,"q",q,"e",e)
[x,y,s,z,info,status] =ecos(c,G,h,dims,A,b);
//  !--error 10000 
// ecos: Expected type ["constant"] for input argument l at input #4, but got "string" instead.
// at line      56 of function Checktype called by :  
// at line     165 of function ecos called by :  
// [x,y,s,z,info,status] =ecos(c,G,h,dims,A,b);
// at line      17 of exec file called by :    
// exec ecos_l1.sce
