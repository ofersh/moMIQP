/*********************************************
 * OPL 12.8.0.0 Model
 * Author: ofers
 * Creation Date: Jan 16, 2023 at 14:02:19
 *********************************************/
int NumOfObj=...;
range rObjectives=1..NumOfObj;
int MaxNumOfPoints=...;
int n=...;
float w[rObjectives]=...;
float solution[i in 0..MaxNumOfPoints][j in rObjectives]=0;
float epsilon=...;
float C=...;
float c1[1..2*n]=...;
float c2[1..2*n]=...;
float Hess[1..2*n][1..2*n]=...;
int LB=...;
int UB=...;
int run_time = ...;
execute {
    cplex.tilim = run_time; //computation time limit <====== TIME LIMIT
    cplex.epagap = epsilon;
	cplex.epgap = epsilon;
	cplex.epint = epsilon;
	cplex.epopt = epsilon;
	cplex.eprhs = epsilon;
	cplex.polishafterintsol = 1;
	cplex.threads=8;
	cplex.mipdisplay=4;
}
//
main{
	var ofile=new IloOplOutputFile("EPF.dat");
	var ofile2=new IloOplOutputFile("EPS.dat");
	ofile.writeln("iteration="+0+";");
	ofile.writeln("solution="+thisOplModel.solution+";");
	//ofile.writeln(thisOplModel.solution);
	ofile.close();
	
	var src=new IloOplModelSource("modelR0.mod");
	var def=new IloOplModelDefinition(src);
	var cplex = new IloCplex();
	var Data= new IloOplDataSource("miRotEllipse_1.dat");
	var Data_1=new IloOplDataSource("EPF.dat");
	
	var Opl=new IloOplModel(def,cplex); 
	Opl.addDataSource(Data);
	Opl.addDataSource(Data_1);
	Opl.generate();
	cplex.solve();
	Opl.postProcess();
	
	writeln(cplex.getObjValue());
	writeln(Opl.printSolution());
	ofile2.writeln(Opl.printSolution());
	
	Opl.end();
	Data.end();
	Data_1.end();
	def.end();
	cplex.end();
	src.end();
	
	for(var iteration=1;iteration<=thisOplModel.MaxNumOfPoints;iteration++){
		var src1=new IloOplModelSource("modelR1.mod");
		var def1=new IloOplModelDefinition(src1);
		var cplex1 = new IloCplex();
		var Data1= new IloOplDataSource("miRotEllipse_1.dat");
		var Data11=new IloOplDataSource("EPF.dat");
		
		var Opl1=new IloOplModel(def1,cplex1);
		Opl1.addDataSource(Data1);
		Opl1.addDataSource(Data11);
		Opl1.generate();
		cplex1.solve();
		Opl1.postProcess();
		
		writeln(cplex1.getObjValue());
		writeln(Opl1.printSolution());
		ofile2.writeln(Opl1.printSolution());
		Opl1.end();
		Data1.end();
		Data11.end();
		def1.end();
		cplex1.end();
		src1.end();
	}
	ofile2.close();		
}