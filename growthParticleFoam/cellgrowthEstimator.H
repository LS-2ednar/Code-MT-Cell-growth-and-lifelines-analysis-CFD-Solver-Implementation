Info<< "\n Estimating new Particles to place \n" << endl;

// Calculate when to add new particles & generate ManuelInjection File 

// ------------------------------------------------------------
// LOOK UP THINGS FROM A FILE DICTIONARY AND USING STORED DATA
// ------------------------------------------------------------	

// Accesssing costum dictionary
dictionary CellDict;
const word CellDictName("cellDict");
    
// Input-Output object which holds the path to the dict and its name 
IOobject CellDictIO
(
CellDictName,				// name of the file
mesh.time().constant(), 	// path to wehere the file is located
mesh,						// reference to mesh needed by the constructor
IOobject::MUST_READ			// reading is required
);

// Check for the availability of the dictionary and if it follows the OF format
if (!CellDictIO.typeHeaderOk<dictionary>(true))
     FatalErrorIn(args.executable()) << "specified dictionary cannot be opend " << CellDictName << exit(FatalError);

// Dictionary object initialisation
CellDict = IOdictionary(CellDictIO);

// Read information form main part of dictionary using standard C++ stringstream syntax
int InitParticles;
CellDict.lookup("InitParticles") >> InitParticles;
	
float Substrate;
CellDict.lookup("Substrate") >> Substrate;
    
float MuMax;
CellDict.lookup("MuMax") >> MuMax;
	
float KS;
CellDict.lookup("KS") >> KS;
	
float qs;
CellDict.lookup("qs") >> qs;
	

// used variables for CellDensity Calculation
double exp(double x);	// initializing the exponential function

// -------------------------------------------------------
// GENERATING AN OUTPUTFILE AND FOLDER TO WORK WITH LATER
// -------------------------------------------------------
    
// Creating custom directory and wirte output file
fileName outputDir = mesh.time().path()/"postProcessing";
mkDir(outputDir);

// Outputfile Pointer
autoPtr<OFstream> outputFilePtr;
    
// Openfile in new created directory
outputFilePtr.reset(new OFstream(outputDir/"CellGrowth.csv"));

// Write header to file
outputFilePtr() << "Time,Cells,Substrate" << endl;
    
// initialize MuT and ST for later use in loop
float MuT;
MuT = MuMax ;
float ST;
ST = Substrate;
int PT;
PT = InitParticles;

// Calculation of Cell values and adding time and CellDensity to outputfile
for (float i = 0; i< 7200/0.0002; i++)	
{		
    //Console Output
	float t = i*0.0002;
	
	// Appending Values to csv file
    outputFilePtr() << runTime.timeName() << "," <<  round(InitParticles*exp(MuT*t)) << "," << ST << endl;
	
	//update MuT, PT and, ST
	MuT = MuMax*(ST/(ST+KS));
	PT  = round(InitParticles*exp(MuT*t));
	ST  = Substrate-(PT*qs);
		
}

return 0;
