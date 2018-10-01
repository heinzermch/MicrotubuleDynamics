#include <iostream>
#include <string>
#include <list> // gotta love the lists
#include <ctime> // time, used for random val init
#include <cstdlib>
#include <cmath> // for log
#include <fstream> // file stream
#include <stdio.h> //remove a file
#include <sstream> // string stream, used for conversion from string to int

#define PRINTMT(a) ((a == 1) ? "+" : "-")
#define ISGTP(a) ((a == 1) ? true : false)
#define PRINTEVENTTYPE(a) ((a == 0) ? "attach" : (a==1) ? "detach" : "hydrolize")

using namespace std;

const double LAMBDA = 40.4;  // rate of attachment
const double H = 0.3;         // rate for Hydrolysis
const double K_MINUS_FACTOR = 900; // Also used as variable kSlope
const double K_MINUS_T = 0.0976;       // rate for detachment for GDP
const double K_MINUS_D = K_MINUS_FACTOR * K_MINUS_T;       // rate of detachment for GTP
const double K_PLUS_EQ = 6.0; // K plus equilibrium rate
const double TIME_LIMIT = 100; // stopping time
const int START_MT_SIZE_ARRAY = 20000; // starting size of the MT
const int PRE_CALC_TABLE_MAX_ENTRY = 1000; // With which precision will the CDF be precalculated
const int MAX_GTP = 3; // Max number of GTP in one vertical line
const double PRECISION = 0.001; // Precision of the inverse calculated value
const double START_MESUREMENT_TIME = 1000.01234;
const double PRINT_INTERVAL = 10.0;
const int START_MT_SIZE = 10000;
const int EVENT_HISTORY_ARRAY_SIZE = 10;

enum event_type_t {ATTACH = 0, DEPOLYMERISATION = 1, HYDROLIZE = 2};

typedef double (*cdf)(double);

double lambda, d, kMinusT, kMinusD, time_limit = 0;

int mtArraySize = START_MT_SIZE_ARRAY;        //size of the MT array
int headPos = 0;            //Where do the GTP tubuli start?
int mtSize = 0;             //The total size of the MT
int* mt = new int[START_MT_SIZE_ARRAY]; // Microtubuli array, -1 is value for empty space, otherwise its the number of GTP molecules
int* eventHistory = new int[EVENT_HISTORY_ARRAY_SIZE]; // Keeping track of events to register catastrophies
int eventHistoryPos = 0;
int growthState = 0;
bool inCatastrophe = false;


// Precalculated value tables for the three actions
double precAttT[PRE_CALC_TABLE_MAX_ENTRY+1] = {};
double precDepT[MAX_GTP+1][PRE_CALC_TABLE_MAX_ENTRY+1] = {};
double precHydT[MAX_GTP][PRE_CALC_TABLE_MAX_ENTRY+1] = {};

// A class to store infromation about an event
class MTEvent {
    int pos;
    double time;
    event_type_t type;
    int width, height, gtps;
public:

    MTEvent (event_type_t e, double t, int p, int nb){
        type = e;
        time = t;
        pos = p;
        gtps = nb;
    };

    MTEvent (event_type_t e, double t, int p){
        type = e;
        time = t;
        pos = p;
        gtps = 0;
    };

    event_type_t getType(){
        return type;
    };
    double getTime(){
        return time;
    };
    int getPos(){
        return pos;
    };
    int getGTPS(){
        return gtps;
    };
};

// A class to store information about the catastrophe, gets used at the beginning and at the end of a catastrophie
class CatastropheEvent {
    int start;
    int stop;
    double time;
public:

    CatastropheEvent (int s, double beginT){
        start = s;
        time = beginT;
    };

    void setStop(int currS, double currentTime){
        stop = currS;
        time = currentTime - time;
    };

    int getLength(){
        return start - stop;
    };

    double getDuration(){
        return time;
    };
};

list<CatastropheEvent> catastrophies;


// assume value is between 0 and 1, containig the borders
// The real table size is maxTablePos+1, but the table contains values from 0 to maxTablePos
// means table[maxTablePos] is valid
// consider two points, (x1,y1), (x2,y2) and the poin to be interpolated (x3,y3)
// normalize x1 = 0, x2 = 1;
// y3 = (y2-y1)*x3 + y1
double interpolate(double table[], double value, int maxTablePos){
    int lowerPos = (int) (value*(maxTablePos+0.0));
    int higherPos = lowerPos+1;
    if(higherPos > maxTablePos)
        return table[lowerPos];
    return (table[higherPos]-table[lowerPos])*(value*(maxTablePos+0.0)-lowerPos)+table[lowerPos];
}

double precAtt(double t){
    return interpolate(precAttT, t, PRE_CALC_TABLE_MAX_ENTRY);
}

double precDep(int gtp, double t){
    return interpolate(precDepT[gtp], t, PRE_CALC_TABLE_MAX_ENTRY);
}

double precHyd(int gtp, double t){
    return interpolate(precHydT[gtp-1], t, PRE_CALC_TABLE_MAX_ENTRY);
}


// CDF for the case TTD
double fa(double t, double lt, double ld){
    return -1.0/((lt-ld)*(lt-ld))*(-lt*lt*exp(-ld*t)-ld*exp(-lt*t)*(ld+lt*(ld*t-2)-lt*lt*t));
}

// CDF for the Gamma distribution, case TTT or DDD
double fTTT(double t, double l){
    return (1.0/(2.0))*exp(-l*t)*(l*l*t*t+2*l*t+2);
}

// CDF for the cases TTD, TDT, DTT
double fTTD(double t, double lt, double ld){
    return fa(t,lt,ld);
}

// CDF for the cases DDT,DTD,TDD
double fTDD(double t, double lt, double ld){
    return fa(t,ld,lt);
}


// CDF for the Gamma distributions, case k=1
double CDFG1(double t, double l){
    return exp(-l*t);
}

// CDF for the Gamma distributions, case k=2
double CDFG2(double t, double l){
    return exp(-l*t)*(l*t+1);
}

// CDF for the Gamma distributions, case k=3
double CDFG3(double t, double l){
    return (1.0/2.0)*exp(-l*t)*(l*l*t*t+2*l*t+2);
}

double CDFHydro1(double t){
    return CDFG1(t,H);
}

double CDFHydro2(double t){
    return CDFG2(t,H);
}

double CDFHydro3(double t){
    return CDFG3(t,H);
}

// Returns a function pointer to the appropriate probability distribution for an hydrolizion event
cdf getGammaCallerForHydro(int k){
    switch (k) {
        case 1:
            return CDFHydro1;
        case 2:
            return CDFHydro3;
        case 3:
            return CDFHydro3;
        default:
            return 0;
    }
}

double CDFDepolyTTT(double t){
    return CDFG3(t,kMinusT);
}

double CDFDepolyTTD(double t){
    return fTTD(t, kMinusD, kMinusT);
}

double CDFDepolyTDD(double t){
    return fTDD(t,kMinusD,kMinusT);
}

double CDFDepolyDDD(double t){
    return CDFG3(t,kMinusD);
}

// Returns a function pointer to the appropriate probability distribution for an depolymerisation event
cdf getDepolyCaller(int gtp){
    switch (gtp) {
        case 0:
            return CDFDepolyDDD;
        case 1:
            return CDFDepolyTDD;
        case 2:
            return CDFDepolyTTD;
        case 3:
            return CDFDepolyTTT;
        default:
            return 0;
    }
}


double CDFAttachG3(double t){
    return CDFG3(t,lambda);
}

// Recursive part of the zero search
double inverseValueR(double (*f)(double), double proba, double low, double high){
    double middle = (high-low)/2.0 + low;
    if(abs(f(middle)-proba) < PRECISION)
        return middle;
    if(f(middle)-proba <= 0)
        return inverseValueR(f, proba, low, middle);
    else
        return inverseValueR(f, proba, middle, high);
}

// Start method for the zero search of a function f
double inverseValue(double (*f)(double), double proba){
    double lower, upper, middle;
    lower = 0;
    upper = 1.31415;
    // lower bound sufficienly close
    if(f(lower)-proba < PRECISION)
        return lower;
    // find an upper bound which is below the zero
    while(f(upper)-proba > 0){
        lower = upper;
        upper +=1;
    }
    if(abs(f(upper)-proba) < PRECISION)
        return upper;
    middle = (upper-lower)/2;
    //cout << "lower " << lower << " middle " << middle << " upper " << upper << endl;
    if(f(middle)-proba <= 0)
        return inverseValueR(f, proba, lower, middle);
    else
        return inverseValueR(f, proba, middle, upper);

}

// Print the current microtubuli on the console
void printMT(){
    for(int i =0; i < mtSize; i++)
            cout << "|" << mt[i];
    cout << ">" << endl;
}

// Print information about the event
void printEvent(MTEvent& event){
    if(event.getType() == HYDROLIZE)
        cout << "Event is: " << PRINTEVENTTYPE(event.getType()) << " at time " << event.getTime() << " and position: " << event.getPos() << " GTPS: " << event.getGTPS();
    else
        cout << "Event is: " << PRINTEVENTTYPE(event.getType()) << " at time " << event.getTime() << " and position: " << event.getPos();
    cout << " headPos: " << headPos << endl;
}

// Print a list of events
void printEvents(list<MTEvent>& events){
    for(list<MTEvent>::iterator iter = events.begin(); iter != events.end(); iter++)
            printEvent(*iter);
}

// Return random number between 0 and 1, including the boundaries
double randomNumber(){
    return ((double)rand())/RAND_MAX;
}

// Return random number between 0 and 1, excluding the boundaries
double randomNumberR(){
    double r = randomNumber();
    while(r == 0.0 || r == 1.0)
        r = randomNumber();
    return r;
}

// Exponentially disributed random number between 0 and infinity
double randomNumberLog(){
    return ((-1)*log(randomNumberR()));
}

// Create a list of all possible events at the current state
MTEvent createEvent(){
    double time = 100000;
    double eventTime = 0;
    int gtps = 0;
    MTEvent event(ATTACH,-1,-1);
    // Hydrolysis events
    for(int i = headPos; i < mtSize; i++){
        gtps = mt[i];
        while(gtps > 0){ // Create event for all possible hydrolization events
            eventTime = precHyd(gtps,randomNumber()); // Have to use the distribution
            if(eventTime < time){
                event = MTEvent(HYDROLIZE,eventTime,i,gtps); // Number of events hydrolisis still unclear
                time = eventTime;
            }
            gtps--;
        }
    }
    // Attach event
    eventTime = precAtt(randomNumber());
    if(eventTime < time){
        event = MTEvent(ATTACH, eventTime,mtSize-1);
        time = eventTime;
    }
    // Depolymerization event
    eventTime = precDep(mt[mtSize-1],randomNumber());
    if(eventTime < time)
        event = MTEvent(DEPOLYMERISATION, eventTime, mtSize-1);
    return event;
}

// Copy the old array into the new one
void resizeMTArray(){
    size_t newSize = mtArraySize * 2;
    int* newArr = new int[newSize];
    for(int i = 0; i < mtArraySize; i++)
        newArr[i] = mt[i];
    for(int i = mtArraySize; i < newSize; i++)
        newArr[i] = -1;
    mtArraySize = newSize;
    delete [] mt;
    mt = newArr;
}


// Apply the chosen event to the Microtubuli
void implementEvent(MTEvent& event){
    switch (event.getType()) {
  case ATTACH:
    // Add new GTPs to the front, if the microtubili is already at the size of the array, resize it
    if(mtSize == mtArraySize){
        resizeMTArray();
    }
    // Add the new GTP polymers
    mt[mtSize] = MAX_GTP;
    if(headPos == mtSize-1 && mt[headPos] == 0)
        headPos++;
    mtSize++;
    break;
  case DEPOLYMERISATION:

    mt[event.getPos()] = -1;
    if(headPos == mtSize-1)
        headPos--;
    mtSize--;

    break;
  case HYDROLIZE:
    mt[event.getPos()] -= event.getGTPS();
    if(event.getPos() == headPos && mt[event.getPos()] == 0)
        while(headPos < mtSize-1){
            headPos++;
            if(mt[headPos] > 0)
                break;
        }
    break;
  default:
    cout << "This should never happen! Programming fail" << endl;
  }
}

// Output the graphdata to the file in data.txt the
void printGraphDataToFile(list<double>& yval, list<double>& xval){
    string dataFileName = "/home/grads/tdpg28/graphData.txt";
    remove(dataFileName);
    ofstream myfile;
    myfile.open(dataFileName);
    list<double>::iterator iterx = xval.begin();
    for(list<double>::iterator itery = yval.begin(); itery != yval.end(); itery++){
        myfile << *itery << "\t" << *iterx << "\n";
        iterx++;
    }
    myfile.close();
}

// Give a system command which prints in the python script
void drawWithPython(){
    std::string command = "python plotfile2.py /home/grads/tdpg28/data1.txt -x=1 -y=2";
    system(command.c_str());
}

//parse arguments from command line, returns true if successful
bool parseArguments(int argc, char *argv[], bool doOutput){
    if(doOutput){
        for(int i = 1; i < argc; i++)
            cout << argv[i] << " ";
        cout << endl;
    }
    double kMinusFactor = K_MINUS_FACTOR;
    double kPlusSlope = K_PLUS_EQ;
    if(argc == 3){
        // Error handeling? Suppose default value is 0
        stringstream(argv[1]) >> d;
        kMinusT = 16.0474*d - 7.3282;
        kMinusD = K_MINUS_FACTOR * kMinusT;
        stringstream(argv[2]) >> time_limit;
        lambda = K_PLUS_EQ*MAX_GTP*d;
        if(lambda == 0 || kMinusT == 0 || d == 0)
            return false;
    }
    else if(argc == 6){
        stringstream(argv[1]) >> d;
        stringstream(argv[2]) >> kMinusD;
        stringstream(argv[3]) >> time_limit;
        stringstream(argv[4]) >> kMinusFactor;
        stringstream(argv[5]) >> kPlusSlope;
        kMinusT =  kMinusD / kMinusFactor;
        lambda = kPlusSlope*MAX_GTP*d;
        if(lambda == 0 || kMinusT == 0 || d == 0 || kMinusFactor == 0)
            return false;
    }
    else{
        cout << "Insufficient arguments detected, using hardcoded arguments" << endl;
        lambda = LAMBDA;
        kMinusT = K_MINUS_T;
        kMinusD = K_MINUS_D;
        time_limit = TIME_LIMIT;
    }
    if(doOutput)
        cout << "Arguments are: k+ " << lambda << " k-t "<< kMinusT<< " k-d " << kMinusD << " time limit " << time_limit << endl;
    return true;
}

void fillInputValues(double input[], double minval, double maxval, int steps){
    double stepsize = maxval/steps;
    double val = minval;
    for(int i = 0; i <= steps;i++){
        input[i] = val;
        val += stepsize;
    }
}

void precalculateValues(double input[], double output[], double (*cdf)(double), int tablesize){
    for(int i = 0; i < tablesize; i++)
        output[i] = inverseValue(cdf, input[i]);
}


// Initialise the microtubuli
void initMT(){
    for(int i = 0; i < START_MT_SIZE; i++)
        mt[i] = MAX_GTP;
    for(int i = START_MT_SIZE; i < mtArraySize; i++)
        mt[i] = -1;
    mtSize = START_MT_SIZE;
    headPos = 0;
}


void initDistributions(){
    double input[PRE_CALC_TABLE_MAX_ENTRY+1] = {};
    fillInputValues(input, 0.0, 1.0, PRE_CALC_TABLE_MAX_ENTRY);

    //Precalculate for attachment
    precalculateValues(input, precAttT, CDFAttachG3 , PRE_CALC_TABLE_MAX_ENTRY+1);

    //Precalculate for depolymerisation
    for(int i = 0; i <= MAX_GTP; i++)
        precalculateValues(input, precDepT[i], getDepolyCaller(i), PRE_CALC_TABLE_MAX_ENTRY+1);

    //Precalcualte for Hydrolision
    for(int i =1; i <= MAX_GTP; i++)
        precalculateValues(input, precHydT[i-1], getGammaCallerForHydro(i) , PRE_CALC_TABLE_MAX_ENTRY+1);
}

// This method detects catastrophies and decides in which state the microtubuli currently is
void detectCatastrophies(MTEvent event, double currentTime){
    if(event.getType() == HYDROLIZE)
        return;
    //Adapt current state and put the current value in the history
    growthState -= eventHistory[eventHistoryPos];
    if(event.getType() == DEPOLYMERISATION)
        eventHistory[eventHistoryPos] = -1;
    else
        eventHistory[eventHistoryPos] = +1;
    growthState += eventHistory[eventHistoryPos];
    eventHistoryPos = (eventHistoryPos+1)%10;
    //See if we are in a catastrophie, if we just entered in one or exited
    if(inCatastrophe && growthState >= 0){ //Exit of a catastrophe
        inCatastrophe = false;
        int stopsAt = 0;
         while(eventHistory[(eventHistoryPos - stopsAt+10)%10] == +1)
            stopsAt = stopsAt -1;
        CatastropheEvent * e = &catastrophies.back();
        e->setStop(mtSize+stopsAt-1, currentTime);
    }
    if(growthState < 0 && !inCatastrophe){ //Entering a catastrophe
        inCatastrophe = true;
        int startsAt = 0;
        while(eventHistory[(eventHistoryPos - startsAt+10)%10] == -1)
            startsAt = startsAt -1;
        CatastropheEvent e(mtSize-startsAt+1, currentTime);
        catastrophies.push_back(e);
    }
}

double runSimulation(bool doGraph, bool doText){
    initMT();             //Initialize the MT
    double currTime = 0;
    list<double> yval;
    list<double> xval;
    double printInt = 0.1;
    double nextPrint = 0.0;
    bool startVal = false;
    int mtStartSize = 0;
    double startMesTime = 1000.01;
    //Initialize the graph data
    if(doGraph){
        yval.push_back(mtSize);
        xval.push_back(0);
    }


    //Going trough the MonteCarlo loop
    while(currTime < time_limit){
        // create and implement the event
        MTEvent event = createEvent();
        implementEvent(event);

        currTime += event.getTime();
        //detectCatastrophies(event, currTime);
        //
        if(currTime > startMesTime && !startVal){
            startVal = true;
            mtStartSize = mtSize;
        }

        // Add the graph info
        if(currTime > nextPrint && doGraph){
            yval.push_back(mtSize);
            xval.push_back(currTime);
            nextPrint += printInt;
        }
        // abort,if the microtubuli has disintegrated
        if(mtSize==0){
            if(doGraph){
                yval.push_back(mtSize);
                xval.push_back(currTime);
            }
            break;
        }
    }
    //If the process had ended in a catastrophe, we need to end the catastrophie
    if(inCatastrophe){
        CatastropheEvent * e = &catastrophies.back();
        e->setStop(mtSize, currTime);
    }
    //Print the graph information to a file and run a python script to display the data
    if(doGraph){
        printGraphDataToFile(xval,yval);
        drawWithPython();
    }
    if(doText){
        cout << "Process finished at " << currTime << " with an MT Size of " << mtSize << endl;
        cout << "Average growth: " << (0.0+mtSize-mtStartSize)/(currTime - startMesTime) << endl;
    }
    return (0.0+mtSize-mtStartSize)/(currTime - startMesTime);
}

// Starts the simulation
int main(int argc, char *argv[])
{
    if(!parseArguments(argc, argv, false))
        return 1;
    //seed random number generator
    srand((unsigned)time(0));
    // Precalculate the distributions
    initDistributions();
    // Run the simulation and output the growth rate
    cout << runSimulation(true,false)*60*0.008 << endl;
    return 0;
}


