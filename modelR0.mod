/*********************************************
 * OPL 12.8.0.0 Model
 * Author: ofers
 * Creation Date: Nov 16, 2022 at 14:02:19
 *********************************************/
int NumOfObj=...;
range rObjectives=1..NumOfObj;
int MaxNumOfPoints=...;
int n=...;
range Dimension=1..n;
float w[rObjectives]=...;
int LB=...;
int UB=...;
dvar float x[d in Dimension] in LB..UB;
dvar int z[d in Dimension] in LB..UB; // <=== INT
dvar float obj[r in rObjectives];
float solution[0..MaxNumOfPoints][rObjectives]=...;
int iteration=...;
int run_time = ...;
float C=...;
float epsilon=...;
float c1[1..2*n]=...;
float c2[1..2*n]=...;
float Hess[1..2*n][1..2*n]=...;
execute {
    cplex.tilim = run_time; //computation time limit
    cplex.epagap = epsilon;
	cplex.epgap = epsilon;
	cplex.epint = epsilon;
	cplex.epopt = epsilon;
	cplex.eprhs = epsilon;
	cplex.polishafterintsol = 1;
	cplex.threads=8;
	cplex.mipdisplay=4;
}    
/******************************************************************************
 * MODEL
 ******************************************************************************/
dexpr float y[i in 1..2*n] = (i in Dimension)?x[i]:z[i-n];
minimize sum(m in rObjectives)(w[m]*obj[m]);
subject to{
	obj[1]==(sum(d1,d2 in 1..2*n) Hess[d1][d2]*(y[d1]-c1[d1])*(y[d2]-c1[d2]))/C; 
	obj[2]==(sum(d1,d2 in 1..2*n) Hess[d1][d2]*(y[d1]-c2[d1])*(y[d2]-c2[d2]))/C; 
}
execute{
	for(var i in rObjectives){
		solution[0][i]=obj[i];
	}
	var ofile=new IloOplOutputFile("EPF.dat");
	ofile.writeln("iteration="+0+";");
	ofile.writeln("solution="+solution+";");
	ofile.close();
}