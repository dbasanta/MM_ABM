package Bone.boneRemodeling_2022May17;

import HAL.GridsAndAgents.*;
import HAL.Gui.*;
import HAL.Interfaces.SerializableModel;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;


import static Bone.boneRemodeling_2022May17.BoneGrid_2022May17.*;
import static HAL.Util.*;
//import static Framework.Util.RGB;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     GRID CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public class BoneGrid_2022May17 extends AgentGrid2D<BoneCell_2022May17> implements SerializableModel {

    ///////////////
    //GRID FIELDS//
    ///////////////

    //Variables to switch on/off treatment
    public static boolean TGFB_INHIBITOR = false;
    public static boolean RANKL_INHIBITOR = false;
    public static boolean BISPHOSPHONATE = false;
    public static boolean BORTEZOMIB = false;
    public static boolean INDIRECT_EFFECT = true;
    public static boolean MYELOMA = false;
    public static boolean EMDR = false;
    public static boolean TREATMENT_ON = false; //this is to control treatment on/off timer in MAIN

    //CLUSTER
    public static boolean PARAM_SWEEP = false; //true; //true; //use when importing parameters to loop through
    public static boolean HEADLESS = false; //true; //use true with cluster


    public final static int BONE = Util.RGB256(255,255,250), MSC = Util.RGB256(135,206,250),
            pOB = Util.RGB256(100,149,237), aOB = Util.BLUE, pOC = Util.RGB256(230,100,130),
            aOC = Util.RED, LINING = Util.RGB256(64,106,151), MM = Util.RGB256(0,128,0) ;

    //SETUP
    static double MinToHour = 60.0;
    public final static double SPACESTEP = 10.0;//um
    public static double TIMESTEP_AGENT = 6.0/MinToHour; //0.1;//hr; //6.0 min; 6.0/60.0 hour
    public final static double N_TIMESTEP_PDE = 60.0*(MinToHour*TIMESTEP_AGENT);//360.0; //Number of diffusion timesteps; 1 dts = 1 sec; 360 dts = 1 ts = 6 min

    //BDF AND RANKL
    double RANKL_productionRate = 2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//2.04e-9//dts; RANKL production rate
    double RANKL_decayRate = -0.35*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //dts; BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE
    double TGFB_productionRate = 2.04e-9*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE;//2.04e-9//dts; RANKL production rate
    static double TGFB_basalRate = 2.04e-11*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //2.04e-11//dts; basal TGFB production rate (not by OC)
    static double TGFB_decayRate = -0.35*(MinToHour*TIMESTEP_AGENT)/N_TIMESTEP_PDE; //dts; BE CAREFUL THIS ISN'T TOO BIG OR ELSE CONCENTRATION GOES NEGATIVE

    //DiffCoef MUST <0.25 for FTCS scheme!
    double RANKL_DiffCoef = 780.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE); //dts
    double TGFB_DiffCoef = 780.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP*N_TIMESTEP_PDE); //dts
    double Extra_TGFBtime = 4320.0/(MinToHour*TIMESTEP_AGENT); //4320 min = 720 ts = 3 days
    static double maxTGFB = 8.7e-10; //1.12e-9; //1.78e-9;//1.21e-9; //6.4e-10
    double maxRANKL = 1.7e-9; //1.02e-9;//1.21e-9;//1.97e-9; //5.8e-10
    double Ts = TGFB_basalRate/(Math.abs(TGFB_decayRate)*maxTGFB); //Basal TGFB

    //CHEMOTAXIS
    double pOC_DiffCoef = 0.1*3.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//3.0*TIMESTEP_AGENT/(SPACESTEP*SPACESTEP);// //slower in marrow than in water
    double MSC_DiffCoef = 0.01*3.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//3.0*TIMESTEP_AGENT/(SPACESTEP*SPACESTEP);//ts //slower in marrow than in water
    double pOB_DiffCoef = 0.01*3.0*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//3.0*TIMESTEP_AGENT/(SPACESTEP*SPACESTEP);//ts //pOB don't wander away from bone
    double pOC_TaxisCoef = 5.0e10*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//5e10 //ts
    double pOB_TaxisCoef = 5.0e11*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//5e11 //ts
    double MSC_TaxisCoef = 5.0e9*(MinToHour*TIMESTEP_AGENT)/(SPACESTEP*SPACESTEP);//4e11 //ts //currently same as MM_TaxisCoef

    //CELL PARAMETERS
    public int TURNOVER_TIME = (int) (2102400.0/(MinToHour*TIMESTEP_AGENT)); //2102400 min = 350400 ts = 4 years
    boolean MMstart = false; //MM recruitment
    public double MM_radius = 80.0 / SPACESTEP;
    public double pmutate = 0.0;

    //PROBABILISTIC RATES (UNIT TIME = 1 HOUR)
    public double MM_DEATH = 1.0 / 7200 * (MinToHour); //*TIMESTEP_AGENT); //7200 min = 5 days; //5760 = 4 days; 4320 min = 3 days
    public double MM_EMDR_DEATH = 1.0 / 7200 * (MinToHour); //*TIMESTEP_AGENT); //7200 min = 5 days; //5760 = 4 days; 4320 min = 3 days
    public double MM_DEATH_BDF = 1.0 / 72000 * (MinToHour); //*TIMESTEP_AGENT); //72000 min = 50 days; 7920 = 5.5 days
    public double MAX_MSC_DIVISION_RATE = 1.0 / 1440 * (MinToHour); //*G.TIMESTEP_AGENT); //1/(1440 min) = 1/(240 ts) = 1/day
    public double MAX_pOB_DIVISION_RATE = 1.0 / 1440 * (MinToHour); //*G.TIMESTEP_AGENT); //1/(1440 min) = 1/(240 ts) = 1/day
    public double MAX_RESISTANT_DIVISION_RATE = 1.0 / 2880 * (MinToHour); //*G.TIMESTEP_AGENT); //1/(1440 min) = 1/(240 ts) = 1/day
    public double MAX_MM_DIVISION_RATE = 1.0 / 1440 * (MinToHour); //*G.TIMESTEP_AGENT); //1/(1440 min) = 1/(240 ts) = 1/day


    //INHIBITORS AND TREATMENT
    double TGFBi=1.0;//1.5; //fold change to production rate
    double RANKLi=1.0;//1.0; //fold change to production rate
    double dose = 1.0; //1.0 ;//1.0; //1.0
    double MM_DEATH_BTZ_FACTOR = 1.5;
    public int Tx_Interval = (int) (30240.0/(MinToHour*TIMESTEP_AGENT)); //(5760.0/(TIMESTEP_AGENT)); //5760 = 4 day; //30240 min = 21 days
    public int Tx_Duration = (int) (30240.0/(MinToHour*TIMESTEP_AGENT));//(20160.0/(TIMESTEP_AGENT)); //1440 = 1 day; 4320 = 3 days; 525600 = 365 days; 20160 min = 14 days
    public int Tx_Total = (int) (525600.0/(MinToHour*TIMESTEP_AGENT));//(525600.0/(MinToHour*TIMESTEP_AGENT)); // 525600 = 365 days
    public int Start_Time=0;
//    public int Cycles = 0;


    //MODEL TESTS
    public boolean pOBadv = true;//false;
    public boolean BDFadv = true;//false;
    public double aOC_scale = 1.0;
//    public double BTZ_scale = 1.0;
//    public int curI = 1;

    public int BMSCpop;
    public double MarrowArea;
    public int pOCpop;

    double MaxRdiff; //max relative change in RANKL
    double MaxTdiff; //max relative change in TGFB
    double rmax;//=0; //static
    double tmax;//=0;
//    double bmax;//=0;
    double convert_to_days = (MinToHour*TIMESTEP_AGENT)/(60.0*24.0); //1 ts = 6 min = 1/240 day
    int pOB_ID_counter;
    int aOC_ID_counter;
    int R_MM_ID_counter;
    int count_BA = 0;
    int init_BA = 0;
    int last_event;
    int Nts = (int) ((4.0*365.0*24.0*60.0)/(MinToHour*TIMESTEP_AGENT));
    public Rand rn;
    public PDEGrid2D RANKL;
    public PDEGrid2D TGFB;
    public ArrayList<BoneCell_2022May17> MSC_List = new ArrayList<>();
    public ArrayList<BoneCell_2022May17> aOB_List = new ArrayList<>();
//    public ArrayList<BoneCell_2022May17> MM_List = new ArrayList<>();
    public int MM_Division_List_Size = 10;//10 timesteps = 60 min = 1 hr;
    public ArrayList<ArrayList<Integer>> MM_Division_List = new ArrayList<ArrayList<Integer>>(MM_Division_List_Size);
    public ArrayList<Integer> MM_Division_Indices = new ArrayList<Integer>();
    public ArrayList<ArrayList<Integer>> MM_Death_List = new ArrayList<ArrayList<Integer>>(MM_Division_List_Size);
    public ArrayList<Integer> MM_Death_Indices = new ArrayList<Integer>();


    public int[] moveHood = VonNeumannHood(true); //Have option of no movement
    public int[] VNHood = VonNeumannHood(false); //4 neighbors
    public int[] VNHood2 = VonNeumannHood(false); //4 neighbors (use if need VNHood twice in loop)
    public int[] MHood = MooreHood(false); //8 neighbors

    public ArrayList<Integer> InitBoneList = new ArrayList<>();
    public ArrayList<BoneCell_2022May17> AllBoneList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<BoneCell_2022May17> LiningList = new ArrayList<>(); //This list is used to randomly determine where remodeling event occurs
    public ArrayList<BoneCell_2022May17> initEventList = new ArrayList<>(); //This list stores which cells express RANKL when event occurs
    public ArrayList<BoneCell_2022May17> tempEventList = new ArrayList<>(); //This list temporarily stores the nearest bone-lining cell neighbors
    public ArrayList<BoneCell_2022May17> recursionList = new ArrayList<>(); //List of LINING cells that we loop through to find nearest bone-lining cell neighbors
    public ArrayList<Integer> tempLiningList = new ArrayList<>();
    public ArrayList<Integer> tempOBList = new ArrayList<>();

    double[] RANKLvals = new double[xDim*yDim];//new double[xDim*yDim];
    double[] TGFBvals = new double[xDim*yDim];//new double[xDim*yDim];
    List<Object> BA_index = new ArrayList<>();
    boolean [][] exposedBone = new boolean[xDim][yDim];
    boolean [][] extraTGFB = new boolean[xDim][yDim];
    double [][] TGFBtimer = new double[xDim][yDim];
    int [][] aOC_depth = new int[xDim][yDim];
//    int [][] BMU_ID = new int[xDim][yDim];
    int [][] resorbedBone = new int[xDim][yDim];


    //public int eOpt;
    FileIO out;
    FileIO RANKLout;
    FileIO TGFBout;
    FileIO BAout;
    FileIO OC_TGFBout;
    FileIO pOBborn;
    FileIO pOBdiff;
    FileIO aOBdeath;
    FileIO MM_birth_death;
    FileIO R_MM_clone;
    FileIO clones;
    FileIO InitialBone;
//    FileIO params;


    //This is important for serializable model
    @Override
    public void SetupConstructors(){
        this._PassAgentConstructor(BoneCell_2022May17.class);
    }

    ////////////////////
    //GRID CONSTRUCTOR//
    ////////////////////

    public BoneGrid_2022May17(int xDim, int yDim, Rand rn, String Bone_FileName) {
        super(xDim, yDim, BoneCell_2022May17.class,true,true);
        this.rn = rn;

        //Create 2D PDE Grid for RANKL and boundary condition
        RANKL = new PDEGrid2D(xDim, yDim,true,true); //This assumes PERIODIC BOUNDARY CONDITION (wrapX=TRUE,wrapY=TRUE)
        TGFB = new PDEGrid2D(xDim, yDim,true,true); //This assumes PERIODIC BOUNDARY CONDITION (wrapX=TRUE,wrapY=TRUE)

        //Create file to record output
//        out = new FileIO(outFileName, "w");
//        out.Write("Timestep" + "," + "BONE" + "," + "pOB" + "," + "aOB" + "," + "pOC" + "," + "aOC" + "," + "MSC" + "," + "LINING" + "," + "MM" + "\n");
//
//        RANKLout=new FileIO(RANKL_FileName,"w");
//        TGFBout=new FileIO(TGFB_FileName,"w");
//        BAout=new FileIO(BA_FileName,"w");
//        OC_TGFBout=new FileIO(OC_TGFB_FileName,"w");
        //OC_TGFBout.Write("Timestep" + "," + "cell_ID" + "," + "N_tick" + "," + "TGFB" + "\n");

        InitialBone=new FileIO(Bone_FileName, "r");

    }

    /////////////////////////
    //GRID METHODS///////////
    /////////////////////////
    //1. InitBone          //
    //2. RemodelingEvent   //
    //3. InitRANKL         //
    //4. ModelStep         //
    //5. CollectLINING     //
    //6. Draw              //
    //7. DrawRANKL         //
    //8. RecordRANKL       //
    //9. RecordOut         //
    /////////////////////////


    //sample from a bounded  distribution
    public double boundedGaussian(double mean, double dev, double min, double max) {
        double gauss = rn.Gaussian(0, 1);
        double val = dev * gauss + mean;
        while (val > max || val < min) {
            gauss = rn.Gaussian(0, 1);
            val = dev * gauss + mean;
        }
        return val;
    }

    public void newFileIO (String projPath, String mode) {

        pOBborn = new FileIO(projPath + "pOBborn.csv", mode);
        pOBdiff = new FileIO(projPath + "pOBdiff.csv", mode);
        aOBdeath = new FileIO(projPath + "aOBdeath.csv", mode);
        MM_birth_death = new FileIO(projPath + "MM_birth_death.csv", mode);
        R_MM_clone = new FileIO(projPath + "R_MM_clone.csv", mode);
        clones = new FileIO(projPath + "clones.csv", mode);
        out = new FileIO(projPath + "PopOut.csv", mode);
        RANKLout = new FileIO(projPath + "RANKLout.csv", mode);
        TGFBout = new FileIO(projPath + "TGFBout.csv", mode);
        BAout = new FileIO(projPath + "BAout.csv", mode);
        OC_TGFBout = new FileIO(projPath + "OC_TGFBout.csv",mode);
//        params = new FileIO(projPath + "params.csv",mode);


        if(mode=="w") {
            pOBborn.Write("cell type" + "," + "cell_ID" + "," + "timestep born" + "," + "RANKL" + "," + "TGFB" + "\n");
            pOBdiff.Write("pOB_ID" + "," + "timestep diff" + "," + "pOB age" + "," + "TGFB" + "\n");
            aOBdeath.Write("time" + "," + "cell type" + "," + "lifespan" + "," + "cell_ID" + "," + "N_tick" + "," + "count_BA" + "," + "TotResorb" + "\n");
            MM_birth_death.Write("action" + "," +  "timestep" +  "," + "nearest bone" + "," + "RANKL" + "," + "TGFB" + "\n");
            R_MM_clone.Write("Timestep" +  "," + "TREATMENT_ON" + "," + "dose" + "," + "TGFBthresh" + "," + "MSC_pOB" + "\n");
            clones.Write("Timestep" + "\n");
            out.Write("Timestep" + "," + "BONE" + "," + "pOB" + "," + "aOB" + "," + "pOC" + "," + "aOC" + "," + "MSC" + "," + "LINING" + "," + "S_MM" + "," + "R_MM" + "," + "TREATMENT_ON" + "," + "BORTEZOMIB" + "," + "MYELOMA" + "\n");
//            params.Write("TGFB_INHIBITOR" + "," + "RANKL_INHIBITOR" + "," + "BISPHOSPHONATE" + "," + "MYELOMA" + "," + "TGFBi" + "," + "RANKLi" + "\n");
//            params.Write(TGFB_INHIBITOR + "," + RANKL_INHIBITOR + "," + BISPHOSPHONATE + "," + MYELOMA + "," + TGFBi + "," + RANKLi + "\n");
        }

    }

    public void closeFileIO () {
        pOBborn.Close();
        pOBdiff.Close();
        aOBdeath.Close();
        MM_birth_death.Close();
        R_MM_clone.Close();
        clones.Close();
        out.Close();
        RANKLout.Close();
        TGFBout.Close();
        BAout.Close();
        OC_TGFBout.Close();
//        params.Close();
    }

    public void SetParams(int prow, ArrayList<String> param_list){
        //returns an array list of all lines from the file as strings

        String[] split_param_list = param_list.get(prow).split(",");

            //RANKL_INHIBITOR=Boolean.parseBoolean(split_param_list[0]);
            //RANKLi=Double.parseDouble(split_param_list[1]);
            //BISPHOSPHONATE = Boolean.parseBoolean(split_param_list[0]);
            //BORTEZOMIB = Boolean.parseBoolean(split_param_list[0]);
            //dose = Double.parseDouble(split_param_list[1]);
            // MM_DEATH = 1.0 / Double.parseDouble(split_param_list[1]) * TIMESTEP_AGENT;
            //MM_radius = Double.parseDouble(split_param_list[1]) / SPACESTEP;
            //pOBadv = Boolean.parseBoolean(split_param_list[1]);
            //aOC_scale = Double.parseDouble(split_param_list[1]);
            //BDFadv = Boolean.parseBoolean(split_param_list[1]);
            //MM_Division_List_Size = Integer.parseInt(split_param_list[1]);
            //BORTEZOMIB = Boolean.parseBoolean(split_param_list[1]);
            //dose = Double.parseDouble(split_param_list[2]);
            //MAX_RESISTANT_DIVISION_RATE = 1.0 / Double.parseDouble(split_param_list[5]) * (MinToHour);
            //MM_DEATH_BDF = 1.0 / Double.parseDouble(split_param_list[5]) * (MinToHour);
            //MM_DEATH_BTZ_FACTOR = Double.parseDouble(split_param_list[5]);
            //MAX_MM_DIVISION_RATE = 1.0 / Double.parseDouble(split_param_list[5]) * (MinToHour);
            //MM_EMDR_DEATH= 1.0 / Double.parseDouble(split_param_list[5]) * (MinToHour);
            //pmutate = Double.parseDouble(split_param_list[2]);

            MYELOMA = Boolean.parseBoolean(split_param_list[0]);
            BORTEZOMIB = Boolean.parseBoolean(split_param_list[1]);
            pmutate = Double.parseDouble(split_param_list[2]);
            EMDR = Boolean.parseBoolean(split_param_list[3]);
            dose = Double.parseDouble(split_param_list[4]);
//            Tx_Duration = (int) (Double.parseDouble(split_param_list[5])/(MinToHour*TIMESTEP_AGENT));

            //MM_DEATH_BDF = 1.0 / Double.parseDouble(split_param_list[5]) * TIMESTEP_AGENT;
            //INDIRECT_EFFECT = Boolean.parseBoolean(split_param_list[4]);
            //TIMESTEP_AGENT = Double.parseDouble(split_param_list[5])/MinToHour; //this doesn't work because other parameters depend on it

    }

    public void InitBone() {

//  FOR RECTANGULAR BONE
//      int xDimBone=60; //px
//      int yDimBone=50;//px
//      int TrSp=100;//Trabecular Spacing (px)
//      int xstart = TrSp/2;//xDim / 5; //initial bone placement
//      int xend = xstart+xDimBone;//4 * xDim / 5; //initial bone placement
//      int ystart = TrSp/2;//2 * yDim / 5; //initial bone placement
//      int yend = ystart+yDimBone;//3 * yDim / 5; //initial bone placement

        //Place bone
        //Bone is placed in center of grid so that 12% of total area is bone.
//        for (int xi = xstart; xi < xend; xi++) {
//            for (int yi = ystart; yi < yend; yi++) {
//                if (xi == xstart || xi == xend - 1 || yi == ystart || yi == yend - 1) {
//                    NewAgentSQ(xi, yi).type = LINING;
//                    GetAgent(xi,yi).Init();
//                    GetAgent(xi,yi).liningAge = TURNOVER_TIME;
//                    LiningList.add(GetAgent(xi, yi));
//                    GetAgent(xi,yi).liningAge = TURNOVER_TIME;
//                } else {
//                    NewAgentSQ(xi, yi).type = BONE;
//                    GetAgent(xi,yi).Init();
//                }
//
//                InitBoneList.add(GetAgent(xi,yi).Isq());
//            }
////            init_BA=InitBoneList.size();
//        }
//        init_BA=InitBoneList.size();

//  FOR IRREGULAR BONE
        int xinit, yinit;
        ArrayList<String> input_data = InitialBone.Read();
        String[] split_input_data =input_data.get(0).split(",");

        //Place bone
        for (int index=1; index<split_input_data.length; index++){
            NewAgentSQ(Integer.parseInt(split_input_data[index])).type=BONE;
            GetAgent(Integer.parseInt(split_input_data[index])).Init();
            InitBoneList.add(Integer.parseInt(split_input_data[index]));
            AllBoneList.add(GetAgent(Integer.parseInt(split_input_data[index])));
        }
        for (int index=1; index<split_input_data.length; index++){
            if(GetAgent(Integer.parseInt(split_input_data[index])).MarrowInHood()==true){
                GetAgent(Integer.parseInt(split_input_data[index])).type=LINING;
                GetAgent(Integer.parseInt(split_input_data[index])).liningAge = TURNOVER_TIME;
                LiningList.add(GetAgent(Integer.parseInt(split_input_data[index])));
            }
        }
        init_BA=InitBoneList.size();
        MarrowArea = (xDim*yDim)-init_BA;//(xDimBone*yDimBone); //0.12 Bone, 0.88 Marrow
        pOCpop = (int) (0.028 * MarrowArea); //0.028 NEED TO PARAMETERIZE!
        BMSCpop = (int) (0.0001 * MarrowArea); //BMSC consist of 0.001%-0.01% of bone marrow

        //Place pOC; random distribution
        for (int i = 0; i < pOCpop; i++) {
            do {
                xinit = rn.Int(xDim);
                yinit = rn.Int(yDim);
            }
            while (PopAt(xinit, yinit) > 0);
            //Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
            NewAgentSQ(xinit, yinit).type = pOC; //Initial stroma
        }

        //Place BMSC; random distribution
        for (int i = 0; i < BMSCpop; i++) {
            do {
                xinit = rn.Int(xDim);
                yinit = rn.Int(yDim);
            }
            while (PopAt(xinit, yinit) > 0);
            //Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
            NewAgentSQ(xinit, yinit).type = MSC; //Initial stroma
            MSC_List.add(GetAgent(xinit,yinit));
        }
    }

    public void RemodelingEvent(ArrayList<BoneCell_2022May17> initEventList, int time){
        //BoneCell_2022May17 initEvent = LiningList.get(rn.Int(LiningList.size())); //todo: uncomment this when want stochastic; only need to use for single event case
        //BoneCell_2022May17 initEvent = GetAgent(399,149);//TODO: UNCOMMENT THIS WHEN DON'T WANT STOCHASTIC
        //initEvent.RANKL_on=true;
        //initEventList.add(initEvent);

        for(int i=0;i<initEventList.size();i++){
            GetAgent(initEventList.get(i).Isq()).type = BONE;//TODO: may want to dispose instead of change to bone to preserve bone density
            GetAgent(initEventList.get(i).Isq()).RANKL_on=true;
            GetAgent(initEventList.get(i).Isq()).eventtime=time;
            //LiningList.remove(GetAgent(initEventList.get(i).Isq())); //commented out 8/29
        }
        initEventList.clear();
    }

    public void EndRemodelingEvent(ArrayList<BoneCell_2022May17> initEventList, int time){
        for(int i=0;i<initEventList.size();i++){
                GetAgent(initEventList.get(i).Isq()).type = LINING;//TODO: may want to dispose instead of change to bone to preserve bone density
                GetAgent(initEventList.get(i).Isq()).RANKL_on = false;
                GetAgent(initEventList.get(i).Isq()).MM_RANKL_on = false;
                GetAgent(initEventList.get(i).Isq()).RANKLtimer = 0;
                GetAgent(initEventList.get(i).Isq()).eventtime = time;
        }
        initEventList.clear();
    }


    public ArrayList<BoneCell_2022May17> InitRANKL(ArrayList<BoneCell_2022May17> recursionList) {
        int[] RANKLHood = MooreHood(false);//Needs to be MooreHood so pOC form on staircase
        int[] MarrowHood = VonNeumannHood(false); //To mirror MarrowInHood
        int MM_count_new=0;

        ArrayList<Integer> MarrowList = new ArrayList<>();

        do{
            int rSize = recursionList.size();
            for (int i = 0; i < rSize; i++) { //Number of cells to check Moore neighborhood
                if(initEventList.size()!=5 && !initEventList.contains(recursionList.get(i))){
                    initEventList.add(recursionList.get(i));
                    int eOpt = recursionList.get(i).MapOccupiedHood(RANKLHood); //Number of cells in Moore neighborhood
                    for (int j = 0; j < eOpt; j++) { //Determine if cells are Bone lining and are not already in initEventList
                        if (GetAgent(RANKLHood[j]).type == LINING && GetAgent(RANKLHood[j]).BuriedInBone()==false && !initEventList.contains(GetAgent(RANKLHood[j]))) { //&& GetAgent(RANKLHood[j]).MarrowInHood()==true
                            int mOpt = GetAgent(RANKLHood[j]).MapHood(MarrowHood);
                            int count_new=0;
                            MM_count_new=0;
                            for(int k=0; k<mOpt; k++){
                                if(GetAgent(MarrowHood[k])==null && !MarrowList.contains(MarrowHood[k])) {
                                    MarrowList.add(MarrowHood[k]);
                                    count_new++;
                                } else if(GetAgent(MarrowHood[k])!=null && GetAgent(MarrowHood[k]).type==MM){
                                    MM_count_new++;
                                }
                            }

                            if(count_new>0 || MM_count_new>0) {
                                tempEventList.add(GetAgent(RANKLHood[j]));
                            }
                        }
                    }
                }
            }
            recursionList.clear();
            recursionList.addAll(tempEventList);
            tempEventList.clear();
        } while(initEventList.size()!=5 && recursionList.size()!=0);
        if (initEventList.size() != 5 || (MarrowList.size() < 5 && MM_count_new==0)) { //if there is not enough marrow for aOC to form not due to MM taking up space
            initEventList.clear();
        }
        return initEventList;
    }

    public ArrayList<BoneCell_2022May17> EndRANKL(ArrayList<BoneCell_2022May17> recursionList) {
        int[] RANKLHood = MooreHood(false);//Needs to be same as InitRANKL

        do{
            int rSize = recursionList.size();
            for (int i = 0; i < rSize; i++) { //Number of cells to check Moore neighborhood
                if(initEventList.size()!=5 && !initEventList.contains(recursionList.get(i))){
                    initEventList.add(recursionList.get(i));
                    int eOpt = recursionList.get(i).MapOccupiedHood(RANKLHood); //Number of cells in Moore neighborhood
                    for (int j = 0; j < eOpt; j++) { //Determine if cells are Bone lining and are not already in initEventList
                        if (GetAgent(RANKLHood[j]).RANKL_on==true && !initEventList.contains(GetAgent(RANKLHood[j]))) {
                            tempEventList.add(GetAgent(RANKLHood[j]));
                        }
                    }
                }
            }
            recursionList.clear();
            recursionList.addAll(tempEventList);
            tempEventList.clear();
        } while(initEventList.size()!=5 && recursionList.size()!=0);
        return initEventList;
    }

    public void IterateRANKL(Grid2Ddouble RANKL_xDiffArray, Grid2Ddouble RANKL_yDiffArray){
        //Diffusion of RANKL
        RANKL.Diffusion(RANKL_xDiffArray,RANKL_yDiffArray);

        //RANKL.Diffusion(RANKL_DiffCoef); //This assumes PERIODIC OR NO-FLUX BOUNDARY CONDITION (need to specify Dirichlet BC otherwise)

        //Natural Decay of RANKL
        RANKL.MulAll(RANKL_decayRate);

//        if(RANKL_INHIBITOR){
//            RANKLi=100.0; //1.41
//        } else {
//            RANKLi=1.0; //no fold change
//        }
//        RANKL.MulAll(RANKLi*RANKL_decayRate);

        //Production of RANKL
        for(int x = 0; x < RANKL.xDim; x++) {
            for (int y = 0; y < RANKL.yDim; y++) {
                if (GetAgent(x,y)!=null && GetAgent(x, y).RANKL_on == true && GetAgent(x,y).RANKLtimer<GetAgent(x,y).max_RANKL_on) {
                    //RANKL.Add(x,y,RANKL_productionRate/maxRANKL);
                    if(RANKL_INHIBITOR && TREATMENT_ON) {
                        RANKL.Add(x, y, RANKLi * RANKL_productionRate / maxRANKL);//*(1-RANKL.Get(x,y)));
                    } else {
                        RANKL.Add(x, y, RANKL_productionRate / maxRANKL);//*(1-RANKL.Get(x,y)));

                    }
                }
            }
        }

        //Solved using Euler's method--need small timestep so that concentration doesn't go negative
        if (RANKL.GetMin() < 0) {
            System.out.println("WARNING: Negative RANKL concentration");
        }

        MaxRdiff = RANKL.MaxDeltaScaled(1e-18);
        RANKL.Update(); //This step is necessary to update diffusion each time-step

        //Record RANKL output
        //RecordRANKL(RANKLout);
    }

    public void IterateTGFB(Grid2Ddouble TGFB_xDiffArray, Grid2Ddouble TGFB_yDiffArray) {
        //Diffusion of TGFB
        TGFB.Diffusion(TGFB_xDiffArray,TGFB_yDiffArray);

        //TGFB.Diffusion(TGFB_DiffCoef); //This assumes PERIODIC OR NO-FLUX BOUNDARY CONDITION (need to specify Dirichlet BC otherwise)

        //Natural Decay of TGFB
        TGFB.MulAll(TGFB_decayRate);

        //TGFB.MulAll(TGFB_decayRate);
//        if(TGFB_INHIBITOR && TREATMENT_ON){
//            TGFBi=2.0;//1.41;
//        } else {
//            TGFBi=1.0; //no fold change
//        }
//        TGFB.MulAll(TGFBi*TGFB_decayRate);

        //Production of TGFB
        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if(TGFB_INHIBITOR && TREATMENT_ON){
                    TGFB.Add(x, y, TGFBi * TGFB_basalRate / maxTGFB);
                    if (GetAgent(x, y) != null && GetAgent(x, y).type == aOC && GetAgent(x, y).TGFB_on == true) {
                        TGFB.Add(x, y, TGFBi * TGFB_productionRate / maxTGFB);
                    } else if (GetAgent(x, y) != null && GetAgent(x, y).type == aOC && GetAgent(x, y).TGFB_on == false) {
                        TGFB.Add(x, y, 0.1 * TGFBi * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    } else if (extraTGFB[x][y] == true && TGFBtimer[x][y] < Extra_TGFBtime && GetAgent(x, y) == null) {
                        TGFB.Add(x, y, 0.1 * TGFBi * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    } else if (extraTGFB[x][y] == true && TGFBtimer[x][y] < Extra_TGFBtime && GetAgent(x, y).type == pOB) {
                        TGFB.Add(x, y, 0.01 * TGFBi * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    }
                } else {
                    TGFB.Add(x, y, TGFB_basalRate / maxTGFB);
                    if (GetAgent(x, y) != null && GetAgent(x, y).type == aOC && GetAgent(x, y).TGFB_on == true) {
                        TGFB.Add(x, y, TGFB_productionRate / maxTGFB);
                    } else if (GetAgent(x, y) != null && GetAgent(x, y).type == aOC && GetAgent(x, y).TGFB_on == false) {
                        TGFB.Add(x, y, 0.1 * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    } else if (extraTGFB[x][y] == true && TGFBtimer[x][y] < Extra_TGFBtime && GetAgent(x, y) == null) {
                        TGFB.Add(x, y, 0.1 * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    } else if (extraTGFB[x][y] == true && TGFBtimer[x][y] < Extra_TGFBtime && GetAgent(x, y).type == pOB) {
                        TGFB.Add(x, y, 0.01 * TGFB_productionRate / maxTGFB); //production by mmp less than aOC
                    }
                }
//                } else if (GetAgent(x,y)!=null && GetAgent(x, y).type==MSC){ //BDF from MSC
//                    TGFB.Add(x, y, 0.01*TGFB_productionRate/maxTGFB); //0.5*
//                } else if (GetAgent(x,y)!=null && GetAgent(x, y).type==pOB) { //BDF from pOB
//                    TGFB.Add(x, y, 0.01*TGFB_productionRate/maxTGFB); //0.4*
//                }
            }
        }


        //Solved using Euler's method--need small timestep so that concentration doesn't go negative
        if (TGFB.GetMin() < 0) {
            System.out.println("WARNING: Negative TGFB concentration");
        }

        MaxTdiff = TGFB.MaxDeltaScaled(1e-18);
        TGFB.Update(); //This step is necessary to update diffusion each time-step

        //Record TGFB output
        //RecordTGFB(TGFBout);
    }


    public void GridTick(){
        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if (extraTGFB[x][y] == true ) {
                    TGFBtimer[x][y]++;
                }
                if(TGFBtimer[x][y]>=Extra_TGFBtime && GetAgent(x,y)==null){
                    extraTGFB[x][y]=false;
                    exposedBone[x][y]=false;
                    TGFBtimer[x][y]=0;
                    aOC_depth[x][y] = 0;
                    resorbedBone[x][y] = 0;
                }
            }
        }
    }

    public void ModelStep(FileIO Write_pOBborn,FileIO Write_pOBdiff, FileIO Write_aOBdeath, FileIO Write_MM_birth_death, FileIO Write_R_MM_clone, int time, List<Integer> RMeventTimes, double [] Cell_Counts) {

        //STEP 0: UPDATE GRIDTICK
        GridTick();

        /////////////////////////////////////////////////
        //STEP 1: REACTION-DIFFUSION EQUATION FOR RANKL//
        /////////////////////////////////////////////////
        int i=0;
        double stol = 1.0e-6;//1e-6; //steady-state tolerance

        //RANKL Diffusion coefficient
        Grid2Ddouble RANKL_xDiffArray = new Grid2Ddouble(xDim,yDim);
        Grid2Ddouble RANKL_yDiffArray = new Grid2Ddouble(xDim,yDim);

        RANKL_xDiffArray.SetAll(RANKL_DiffCoef);
        RANKL_yDiffArray.SetAll(RANKL_DiffCoef);

        for(int x=0;x<RANKL.xDim;x++) {
            for (int y = 0; y < RANKL.yDim; y++) {
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING) && GetAgent(x,y).RANKL_on==false) {
                    RANKL_xDiffArray.Set(x, y, 0.1 * RANKL_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (x!=xDim-1 && GetAgent(x+1, y) != null && (GetAgent(x+1, y).type == BONE || GetAgent(x+1, y).type == LINING) && GetAgent(x+1,y).RANKL_on==false){
                    RANKL_xDiffArray.Set(x, y, 0.1 * RANKL_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING) && GetAgent(x,y).RANKL_on==false) {
                    RANKL_yDiffArray.Set(x, y, 0.1 * RANKL_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (y!=yDim-1 && GetAgent(x, y+1) != null && (GetAgent(x, y+1).type == BONE || GetAgent(x, y+1).type == LINING) && GetAgent(x,y+1).RANKL_on==false){
                    RANKL_yDiffArray.Set(x, y, 0.1 * RANKL_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
            }
        }

        //TGFB Diffusion coefficient
        Grid2Ddouble TGFB_xDiffArray = new Grid2Ddouble(xDim,yDim);
        Grid2Ddouble TGFB_yDiffArray = new Grid2Ddouble(xDim,yDim);

        TGFB_xDiffArray.SetAll(TGFB_DiffCoef);
        TGFB_yDiffArray.SetAll(TGFB_DiffCoef);

        for(int x=0;x<TGFB.xDim;x++) {
            for (int y = 0; y < TGFB.yDim; y++) {
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    TGFB_xDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (x!=xDim-1 && GetAgent(x+1, y) != null && (GetAgent(x+1, y).type == BONE || GetAgent(x+1, y).type == LINING)){
                    TGFB_xDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
                if (GetAgent(x, y) != null && (GetAgent(x, y).type == BONE || GetAgent(x, y).type == LINING)) {
                    TGFB_yDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                } else if (y!=yDim-1 && GetAgent(x, y+1) != null && (GetAgent(x, y+1).type == BONE || GetAgent(x, y+1).type == LINING)){
                    TGFB_yDiffArray.Set(x, y, 0.1 * TGFB_DiffCoef); //Originally had RANKL_DiffCoef*0.1
                }
            }
        }


        do {
            IterateRANKL(RANKL_xDiffArray,RANKL_yDiffArray);
            IterateTGFB(TGFB_xDiffArray,TGFB_yDiffArray);
            //Iterate
            i++;
        } while (i<N_TIMESTEP_PDE && Math.max(MaxRdiff,MaxTdiff) >= stol);

        //Record overall max TGFB and RANKL (to know how to normalize)
        if(rmax<RANKL.GetMax()){
            rmax=RANKL.GetMax();
        }
        if(tmax<TGFB.GetMax()){
            tmax=TGFB.GetMax();
        }

        System.out.println("max RANKL "+rmax);
        System.out.println("max TGFB "+tmax);


        /////////////////////////////////
        //STEP 3: ITERATE THROUGH CELLS//
        /////////////////////////////////

        CleanShuffle(rn);
        //ShuffleAgents(rn);
        for (BoneCell_2022May17 c: this) {
            c.CellStep(Write_pOBborn, Write_pOBdiff, Write_aOBdeath, Write_MM_birth_death, Write_R_MM_clone, time, RMeventTimes, Cell_Counts);
        }

        /////////////////////////////////
        //STEP 4: INCREASE AGE OF CELLS//
        /////////////////////////////////
        for (BoneCell_2022May17 c:this){
            c.CellTick();
        }


    }

    ////////////////////////////
    //Full Domain (No zoom-in)//
    ////////////////////////////
//
    public void Draw(UIGrid vis, UILabel days, int i) {
        days.SetText("days: "+i*convert_to_days);
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                if (drawMe != null && drawMe.type==MM && !drawMe.RESISTANT && (drawMe.pOBInCircleHood(drawMe.protect_radius)==true || TGFB.Get(drawMe.Xsq(),drawMe.Ysq())>=drawMe.TGFBthresh)) {
                    vis.SetPix(x,y, Util.GREEN); //MARROW=LIGHT PINK
                } else if (drawMe != null && drawMe.type==MM && drawMe.RESISTANT && (drawMe.pOBInCircleHood(drawMe.protect_radius)==true || TGFB.Get(drawMe.Xsq(),drawMe.Ysq())>=drawMe.TGFBthresh)) {
                    vis.SetPix(x,y, Util.CategorialColor(19)); //MARROW=LIGHT PINK
                } else if (drawMe != null && drawMe.type==MM && drawMe.RESISTANT==true){
                    vis.SetPix(x,y, Util.BLACK);
                } else if (drawMe != null) {
                    //vis.SetPix(x, y, drawMe.color);
                    vis.SetPix(x, y, drawMe.type);
                } else{
                    vis.SetPix(x,y, Util.RGB256(240, 220, 220)); //MARROW=LIGHT PINK
                }
            }
        }
        vis.SetString("Day: "+(int)(i*convert_to_days),1,yDim-1,BLACK,Util.RGB256(240, 220, 220));

    }

    public void DrawDivision(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                vis.SetPix(x, y, Util.BLACK);
                if (drawMe != null) {
//                    if(drawMe.type==MM && drawMe.color!=MM){
//                        vis.SetPix(x, y, drawMe.color);
                    if(drawMe.type==LINING || drawMe.type==BONE) {
                        vis.SetPix(x, y, BONE);//drawMe.type);
                    } else{
                        vis.SetPix(x, y, Util.BLACK);
                    }
                }
            }
        }
        for (int i1 = 0; i1<MM_Division_List.size(); i1++){
            for(int i2=0; i2<MM_Division_List.get(i1).size(); i2++) {
                vis.SetPix(MM_Division_List.get(i1).get(i2), Util.GREEN);
            }
        }
    }

    public void DrawDeath(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                vis.SetPix(x, y, Util.BLACK);
                if (drawMe != null) {
//                    if(drawMe.type==MM && drawMe.color!=MM){
//                        vis.SetPix(x, y, drawMe.color);
                    if(drawMe.type==LINING || drawMe.type==BONE) {
                        vis.SetPix(x, y, BONE);//drawMe.type);
                    } else{
                        vis.SetPix(x, y, Util.BLACK);
                    }
                }
            }
        }
        for (int i1 = 0; i1<MM_Death_List.size(); i1++){
            for(int i2=0; i2<MM_Death_List.get(i1).size(); i2++) {
                vis.SetPix(MM_Death_List.get(i1).get(i2), Util.RED);
            }
        }
    }


    public static int HeatMapLog(double val){
        if ((val>=0) && (val<=0.0001)) {
            return RGB(0,0,0); //black
        } else if ((val>0.0001) && (val<=0.001)) {
            return ColorMap(val,0.0001,0.001,RGB(0,0,0),RGB(0,0,255)); //black->blue
        } else if ((val>0.001) && (val<=0.01)) {
            return ColorMap(val,0.001,0.01,RGB(0,0,255),RGB(0,206,209));//blue->cyan
        } else if ((val>0.01) && (val<=0.1)) {
            return ColorMap(val,0.01,0.1,RGB(0,206,209),RGB(128,0,128));//cyan->purple
        } else if ((val>0.1) && (val<=1)) {
            return ColorMap(val,0.1,1,RGB(128,0,128),RGB(255,0,0)); //purple->red
        } else {
            return RGB(255, 255,  255); //white
        }
    }

    public void DrawRANKL(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else {
                    vis.SetPix(x, y, HeatMapLog(RANKL.Get(x,y)));
                }
            }
        }
    }

    public void DrawTGFB(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                BoneCell_2022May17 drawMe = GetAgent(x,y);
                if(drawMe!=null && drawMe.type==LINING) {
                    vis.SetPix(x, y, BONE);//drawMe.type);
                } else {
                    vis.SetPix(x, y, HeatMapLog(TGFB.Get(x,y)));
                }
            }
        }
    }


    public void RecordRANKL(FileIO PDEwriteHere,int time) {
        int i=0;
        for (int x = 0; x < xDim; x++) { //x<80
            for (int y = 0; y < yDim; y++) { //y<50
                RANKLvals[i]=(RANKL.Get(x, y));//RANKL.Get(x+190,y+125)
                i++;
            }
        }
        PDEwriteHere.Write(time+",");
        PDEwriteHere.WriteDelimit(RANKLvals,",");
        PDEwriteHere.Write("\n");

    }

    public void RecordTGFB(FileIO PDEwriteHere,int time) {
        int i=0;
        for (int x = 0; x < xDim; x++) { //x<80
            for (int y = 0; y < yDim; y++) { //y<50
                TGFBvals[i]=(TGFB.Get(x, y));//RANKL.Get(x+190,y+125)
                i++;
            }
        }
        PDEwriteHere.Write(time+",");
        PDEwriteHere.WriteDelimit(TGFBvals,",");
        PDEwriteHere.Write("\n");

    }

    public void RecordBA(FileIO BAwriteHere,int time) {

        for (BoneCell_2022May17 c : this) {
            if (c.type == BONE || c.type == LINING) {
                BA_index.add(c.Isq());
            }
        }

        BAwriteHere.Write(time+",");
        BAwriteHere.WriteDelimit(BA_index,",");
        BAwriteHere.Write("\n");

    }

    public void RecordClones(FileIO writeHere,int time){
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int[] cts = new int[R_MM_ID_counter];

        for (BoneCell_2022May17 c : this) {
            if(c.type==MM && c.RESISTANT) {
                cts[c.R_MM_ID-1]++;
            }
        }
        //population of one timestep per line
        //writeHere.Write(ct_BONE+","+ct_pOB+","+ct_aOB+","+ct_pOC+","+ct_aOC+","+ct_MSC+","+ct_LINING+"\n");
        writeHere.Write(time+",");
        writeHere.WriteDelimit(cts,",");
        writeHere.Write("\n");
    }

    public void RecordOut(FileIO writeHere,int time, boolean treatment_on, boolean btz, boolean myeloma){
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int[] cts = new int[9];

        for (BoneCell_2022May17 c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            } else if(c.type==pOB){
                //ct_pOB++;
                cts[1]++;
            } else if(c.type==aOB){
                //ct_aOB++;
                cts[2]++;
            } else if(c.type==pOC){
                //ct_pOC++;
                cts[3]++;
            } else if(c.type==aOC){
                //ct_aOC++;
                cts[4]++;
            } else if(c.type==MSC){
                //ct_MSC++;
                cts[5]++;
            } else if(c.type==LINING){
                //ct_LINING++;
                cts[6]++;
            } else if(c.type==MM && !c.RESISTANT){
                cts[7]++;
            } else if(c.type==MM && c.RESISTANT) {
                cts[8]++;
            }
        }
        //population of one timestep per line
        //writeHere.Write(ct_BONE+","+ct_pOB+","+ct_aOB+","+ct_pOC+","+ct_aOC+","+ct_MSC+","+ct_LINING+"\n");
        writeHere.Write(time+",");
        writeHere.WriteDelimit(cts,",");
        writeHere.Write("," + treatment_on + "," + btz + "," + myeloma + "\n");
    }

    public double[] CellCounts(){
        //int BONE = 0, pOB = 1, aOB = 2, pOC = 3, aOC = 4, MSC = 5, LINING = 6, MM = 7;
        double[] cts = new double[9];

        for (BoneCell_2022May17 c : this) {
            if(c.type == BONE){
                //ct_BONE++;
                cts[0]++;
            } else if(c.type==pOB){
                //ct_pOB++;
                cts[1]++;
            } else if(c.type==aOB){
                //ct_aOB++;
                cts[2]++;
            } else if(c.type==pOC){
                //ct_pOC++;
                cts[3]++;
            } else if(c.type==aOC){
                //ct_aOC++;
                cts[4]++;
            } else if(c.type==MSC){
                //ct_MSC++;
                cts[5]++;
            } else if(c.type==LINING){
                //ct_LINING++;
                cts[6]++;
            } else if(c.type==MM && !c.RESISTANT){
                cts[7]++;
            } else if(c.type==MM && c.RESISTANT){
                cts[8]++;
            }
        }
        return cts;
    }

    public boolean RM_CellCounts() {
        //int ct_BONE = 0, ct_pOB = 0, ct_aOB = 0, ct_pOC = 0, ct_aOC = 0, ct_MSC = 0, ct_LINING = 0;
        int cts = 0;
        for (BoneCell_2022May17 c : this) {
            if (c.RANKL_on == true || c.type == aOB || c.type == aOC) {
                cts++;
            }
        }
//        if(cts==0) {
//            for (int x = 0; x < xDim; x++) {
//                for (int y = 0; y < yDim; y++) {
//                    if(exposedBone[x][y]==true){
//                        cts++;
//                    }
//                }
//            }
//        }
        return cts == 0;
    }

    public boolean RM_TGFB_off() {
        int cts=0;
        for (int x = 0; x < xDim; x++) { //x<80
            for (int y = 0; y < yDim; y++) { //y<50
                if (extraTGFB[x][y]==true && TGFBtimer[x][y]<Extra_TGFBtime) {
                    cts++;
                }
            }
        }

        return cts == 0;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                      MAIN                                                      //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static void main(String[] args) {

//        boolean HEADLESS = true; //use true with cluster

        if (HEADLESS) {
            System.setProperty("java.awt.headless", "true");
        }

        int xDim = 160;//500; //px
        int yDim = 150;//250; //px
        String sdf = new SimpleDateFormat("yyyyMMMd").format(new Date());


//        int Nts = 8*365*24*10;//365*24*10; //125 days//30000;//6000;//3*87600; //30000 //Number of timesteps; 1 ts = 6 min, 24000 ts = 100 day

        //int pOB_ID_inc = 0;

        //GridWindow win=new GridWindow(xDim,yDim,3); //Create instance of GridWindow class (imported)
        //GridWindow PDEwin = new GridWindow(xDim,yDim,1);

        UIWindow win = HEADLESS ? null : new UIWindow("Normal Bone Remodeling");
        String fn = "Bone_" + sdf;
        File dir = new File(fn);
        dir.mkdir();

        ///////////////////////
        //FOR PARAMETER SWEEP//
        ///////////////////////

        int param_list_size;
        ArrayList<String> param_list = null;

        if(PARAM_SWEEP) {
            FileIO Params = new FileIO("Bone/params.csv", "r");
            param_list = Params.Read();
            param_list_size = param_list.size();
        }else{
            param_list_size = 2; //this will go through iteration once for pre-defined parameters
        }


        //Define parameters
        for (int prow=1; prow<param_list_size; prow++) {

        for (int sim = 0; sim < 25; sim++) {
            String subfolder = fn + "/Sim" + sim + "_row" + prow + "/";
            dir = new File(subfolder);
            boolean success = dir.mkdir();

            if (success) {


                ////////////////////////////////////////
                //UIGrid                              //
                //compx=number of columns (default 1) //
                //compy=number of rows (default 1)    //
                ////////////////////////////////////////

                //Full domain
                UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
                UIGrid RANKL_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
                UIGrid TGFB_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
                UIGrid MMDiv_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2
                UIGrid MMDeath_vis = new UIGrid(xDim, yDim, 2);//scaleFactor with BDF: 2



                //Sub-domain (zoom-in on remodeling site)
//        UIGrid Cell_vis = new UIGrid(80,50,2,2,1);
//        UIGrid RANKL_vis = new UIGrid(80,50,1);
//        UIGrid TGFB_vis = new UIGrid(80,50,1);

                //Sub-domain (zoom-in on multi-remodeling sites)
//                UIGrid Cell_vis = new UIGrid(350, 100, 2, 2, 1);
//                UIGrid RANKL_vis = new UIGrid(350, 100, 1);
//                UIGrid TGFB_vis = new UIGrid(350, 100, 1);


                UILabel days = new UILabel("days:______________________");


                if (!HEADLESS) {
                    win.AddCol(0, new UILabel("Cells"));
                    win.AddCol(1, days);
                    win.AddCol(0, Cell_vis);
                    win.AddCol(2, new UILabel("RANKL"));
                    win.AddCol(2, RANKL_vis);
                    win.AddCol(2, new UILabel("TGFB"));
                    win.AddCol(2, TGFB_vis);
                    win.AddCol(3, new UILabel("MM Division"));
                    win.AddCol(3, MMDiv_vis);
                    win.AddCol(3, new UILabel("MM Death"));
                    win.AddCol(3, MMDeath_vis);

                    win.RunGui();
                }

                // GIF MAKER

//            String projPath = PWD() + "/Bone/";
                String projPath = subfolder; //PWD() + subfolder;
                GifMaker gm_Cell_vis = new GifMaker(projPath.concat("/").concat("CellVid").concat(".gif"), 100, true);
                GifMaker gm_RANKL_vis = new GifMaker(projPath.concat("/").concat("RANKLVid").concat(".gif"), 100, true);
                GifMaker gm_TGFB_vis = new GifMaker(projPath.concat("/").concat("TGFBVid").concat(".gif"), 100, true);
                GifMaker gm_MMDiv_vis = new GifMaker(projPath.concat("/").concat("MMDivVid").concat(".gif"), 100, true);
                GifMaker gm_MMDeath_vis = new GifMaker(projPath.concat("/").concat("MMDeathVid").concat(".gif"), 100, true);

                // RECORD OUTPUT
//                FileIO pOBborn = new FileIO(projPath + "pOBborn.csv", "w");
//                FileIO pOBdiff = new FileIO(projPath + "pOBdiff.csv", "w");
//                FileIO aOBdeath = new FileIO(projPath + "aOBdeath.csv", "w");
//                pOBborn.Write("cell type" + "," + "cell_ID" + "," + "timestep born" + "," + "RANKL" + "," + "TGFB" + "\n");
//                pOBdiff.Write("pOB_ID" + "," + "timestep diff" + "," + "pOB age" + "," + "TGFB" + "\n");
//                aOBdeath.Write("time" + "," + "cell type" + "," + "lifespan" + "," + "cell_ID"+ "," + "N_tick" + "," + "count_BA" + "," + "TotResorb" + "\n");

                BoneGrid_2022May17 g = new BoneGrid_2022May17(xDim, yDim, new Rand(), "/Users/david/HAL-master/MM/BAout_2020May5_Sim14.csv"); //Create instance of BoneGrid class (grid)
                //Rand(seed:1) when want to reproduce results

                //Set treatment schedule
                ArrayList<Double> Tx_subStart = new ArrayList<>();
                if(g.Tx_Duration==g.Tx_Interval) { //continuous treatment
                    Tx_subStart.add(0.0);//1440.0/TIMESTEP_AGENT); //1 day=1440 minutes
                }else{ //pulsed treatment
                    Tx_subStart.add(0.0);//1440.0/TIMESTEP_AGENT); //1 day=1440 minutes
//                    Tx_subStart.add(3*1440.0/TIMESTEP_AGENT); //4 days
//                    Tx_subStart.add(7*1440.0/TIMESTEP_AGENT); //8 days
//                    Tx_subStart.add(10*1440.0/TIMESTEP_AGENT); //11 days
                }
                double subStart_time=0;


                //Set parameters
                if(PARAM_SWEEP) { //Note: Other parameters dependent on this list will not be updated
                    g.SetParams(prow, param_list);
                }


                //Record Output
                g.newFileIO(projPath, "w");

                //Initialize model
                g.InitBone();

                int NLining = g.LiningList.size();
                int Nevents = g.init_BA / 50;//(14 * 5);//NLining / 5; //on average, OC resorb 14*5 units of bone. This is how many times that happens to resorb init_BA
                int[] RMevents;

//                  RMevents=new int[Nevents];
//                  g.rn.RandomIS(RMevents, 0, g.curI * g.TURNOVER_TIME);
                    RMevents = new int[Nevents * g.Nts / g.TURNOVER_TIME];
                    g.rn.RandomIS(RMevents, 0, g.Nts);

                List<Integer> RMeventTimes = Arrays.stream(RMevents).boxed().collect(Collectors.toList()); //convert array to list
                RMeventTimes.set(0, 0);


                ////////////////
                //Loop through//
                ////////////////

                for (int i = 0; i < g.Nts; i++){//g.Nts
                    System.out.println("Timestep=" + i);

                    if (!HEADLESS) {
                        win.TickPause(100); //slow down visualization
                    }

                    if (i >= ((4.0*365.0*24.0*60.0)/(MinToHour*TIMESTEP_AGENT)) && g.count_BA >= g.init_BA && g.BA_index.isEmpty()) { //TO SAVE BONE
                        g.RecordBA(g.BAout, i);
                    }


//                    if(i<g.TURNOVER_TIME/8.0){
//                        TREATMENT_ON=false;
//                    } else {
//
                    double [] Cell_Counts = g.CellCounts(); //int BONE = 0, pOB = 1, aOB = 2, pOC = 3, aOC = 4, MSC = 5, LINING = 6, MM = 7;

//                    if(TREATMENT_ON == false && MYELOMA == true && (Cell_Counts[7] / (xDim * yDim - (Cell_Counts[0] + Cell_Counts[6])) >= 0.1)){
//                        TREATMENT_ON = true;
//                        g.Start_Time = i;
//                    } else if (TREATMENT_ON == true && MYELOMA == true && i>(g.Start_Time + g.Tx_Duration)){
//                        TREATMENT_ON = false;
//                        g.Nts = i;
//                    }


                if(MYELOMA && BORTEZOMIB) {
                    if (!TREATMENT_ON && ((Cell_Counts[7]+Cell_Counts[8]) / (xDim * yDim - (Cell_Counts[0] + Cell_Counts[6])) >= 0.1) && g.Start_Time == 0) {
                        TREATMENT_ON = true;
                        g.Start_Time = i;
                        subStart_time = i;
                    } else if (g.Start_Time > 0 && ((Cell_Counts[7]+Cell_Counts[8]) / (xDim * yDim - (Cell_Counts[0] + Cell_Counts[6])) >= 0.2)){//i > (g.Start_Time + g.Tx_Total)) {
                        TREATMENT_ON = false;
                        g.Nts = i;
                    } else if (g.Start_Time > 0 && Tx_subStart.contains((double) ((i - g.Start_Time) % g.Tx_Interval))) {
                        subStart_time = i;//(i-g.Start_Time)%g.Tx_Interval;
                        TREATMENT_ON = true;
                    } else if (g.Start_Time > 0 && i > (subStart_time + g.Tx_Duration)) {
                        TREATMENT_ON = false;
                    }
                } else if (!MYELOMA && (BORTEZOMIB || BISPHOSPHONATE || TGFB_INHIBITOR || RANKL_INHIBITOR)){
                    TREATMENT_ON = !(i < g.TURNOVER_TIME / 4.0);
                }


                    g.ModelStep(g.pOBborn, g.pOBdiff, g.aOBdeath, g.MM_birth_death, g.R_MM_clone,i, RMeventTimes,Cell_Counts);

//                 if(g.MM_Death_Indices.size()>0){
//                     System.out.println("Start Here");
//                 }

                    if(i>=g.MM_Division_List_Size) { //TODO: CHECK THIS LOGIC
                        g.MM_Division_List.remove(0);
                        g.MM_Division_List.add(g.MM_Division_Indices);
                        g.MM_Death_List.remove(0);
                        g.MM_Death_List.add(g.MM_Death_Indices);
                    } else {
                        g.MM_Division_List.add(g.MM_Division_Indices);
                        g.MM_Death_List.add(g.MM_Death_Indices);
                    }
                    g.MM_Division_Indices = new ArrayList<Integer>();//.clear();
                    g.MM_Death_Indices = new ArrayList<Integer>();

                    //System.out.println(g.Pop());
                    //System.out.println(g.RANKL.GetMax());
                    //System.out.println(g.RANKL.GetMin());
                    g.IncTick(); //update all cell ages at once
                    g.Draw(Cell_vis, days, i);
                    g.DrawRANKL(RANKL_vis);
                    g.DrawTGFB(TGFB_vis);
                    g.DrawDivision(MMDiv_vis);
                    g.DrawDeath(MMDeath_vis);


//                    if (i%(24.0*60.0/(MinToHour*g.TIMESTEP_AGENT))==0) { //24 ts = 0.1 day
//                        gm_Cell_vis.AddFrame(Cell_vis);
//                        gm_RANKL_vis.AddFrame(RANKL_vis);
//                        gm_TGFB_vis.AddFrame(TGFB_vis);
//                        gm_MMDiv_vis.AddFrame(MMDiv_vis);
//                        gm_MMDeath_vis.AddFrame(MMDeath_vis);
//////                        g.RecordRANKL(g.RANKLout, i);
//////                        g.RecordTGFB(g.TGFBout, i);
//                    }



//            if(i%(Nts/100)==0){
////           if(i==0){
////                Cell_vis.ToPNG("/Users/anna/Desktop/HAL-2019July/Bone"+i+".png");
////                RANKL_vis.ToPNG("/Users/anna/Desktop/HAL-2019July/Bone" + "_RANKL"+i+ ".png");
////                TGFB_vis.ToPNG("/Users/anna/Desktop/HAL-2019July/Bone" + "_TGFB"+i+ ".png");
//            } else if(i%240==0){
//                g.RecordRANKL(g.RANKLout,i);
//                g.RecordTGFB(g.TGFBout,i);
//            }
                    if(i%(24.0*60.0/(MinToHour* TIMESTEP_AGENT))==0) { //240ts = 1 day
                        g.RecordOut(g.out, i,TREATMENT_ON,BORTEZOMIB,MYELOMA);
                        g.RecordClones(g.clones, i);
                    }
                    //Data recording
                    //RANKLout.Write(g.RANKL.GetMin()+"\n");
                }


                //Close FileIO
                g.closeFileIO();

                //Close GifMaker
                gm_Cell_vis.Close();
                gm_RANKL_vis.Close();
                gm_TGFB_vis.Close();
                gm_MMDiv_vis.Close();
                gm_MMDeath_vis.Close();

                //Close UIWindow
                if (!HEADLESS) {
                    win.Close();
                }
                //PDEwin.Close();
            }
        } //END OF MULTIPLE RUNS
        } //end parameter loop
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                     CELL CLASS                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BoneCell_2022May17 extends AgentSQ2Dunstackable<BoneGrid_2022May17> {

    ///////////////
    //CELL FIELDS//
    ///////////////

    public int type;
    public int color;
    int aOC_ID;//defines which aOC unit number
    int pOB_ID; //defines pOB number
    int R_MM_ID; //defines R_MM clone

    boolean TGFB_on = false;
    boolean RANKL_on = false;
    boolean MM_RANKL_on = false;
    int eventtime=(int) ((8.0*365.0*24.0*60.0)/(MinToHour*TIMESTEP_AGENT)); //Nts //defines when cell completed event out of turn so that it is not selected again in same timestep (like the 5 aOC)
    //aOC
    int aOCtick = 0; //counter for aOC to resorb one unit of bone
    int resorbtick = 0; //defines total resorption depth
    int aOCage = 0; //counter for when aOC dies
    boolean aOC_DONE = false;
    //pOB
    int pOBage = 0; //keeps track of pOB age
    int pOBtick = 0; //counter for when pOB differentiates
    //int pOB_ID_label; //need static so it increments every time
    //aOB
    int aOBage = 0; //counter for when aOB dies
    int aOBtick = 0; //counter for aOB to build one unit of bone
    int formtick = 0; //defines total bone formed
    int RANKLtimer = 0; //counter to set max time that RANKL on without event happening
    int NeedBuild=0;
    //Lining
    int liningAge = 0; //counter for when lining old enough for remodeling event
    //MM
    boolean RESISTANT = false;

    public void InitAttributes(int time) {
        color = type;
        TGFB_on = false;
        RANKL_on = false;
        MM_RANKL_on = false;
        eventtime = time; //defines when cell completed event out of turn so that it is not selected again in same timestep (like the 5 aOC)
        //aOC
        aOCtick = 0; //counter for aOC to resorb one unit of bone
        resorbtick = 0; //defines total resorption depth
        aOCage = 0; //counter for when aOC dies
        aOC_DONE = false;
        //pOB
        pOBage = 0; //keeps track of pOB age
        pOBtick = 0; //counter for when pOB differentiates
        //int pOB_ID_label; //need static so it increments every time
        //aOB
        aOBage = 0; //counter for when aOB dies
        aOBtick = 0; //counter for aOB to build one unit of bone
        formtick = 0; //defines total bone formed
        RANKLtimer = 0; //counter to set max time that RANKL on without event happening
        NeedBuild=0;
        //lining
        liningAge = 0;
        //MM
        RESISTANT = false;
    }

    //no longer used
//    int aOC_neighbors;

//    public int[] pOCHood = VonNeumannHood(false);//Check if pOC neighbor is bone or other pOC
//    public int[] boneInHood = VonNeumannHood(false);//check if bone in neighborhood of aOC (should be Von
//                                                                // Neumann because moving diagonal has different propensity)
    //Note: if use Moorehood for boneInhHood it's possible that cells more so that there is no bone in neighborhood
    //so an error is produced. May need to incorporate this scenario into code even for VonNeumannHood.
//    public int[] aOCinHood = MooreHood(false);//check if aOC in neighborhood of bone (can be Moore because
//                                                          //just need to keep aOC connected)

    public int fOpt;
    int pOCx, pOCy;

    //CELL PARAMETERS
    public int aOC_DEATH;// = (int) (20160 / G.TIMESTEP_AGENT); //20160 min = 3360 ts = 14 days
    //    public double MM_DEATH = 1.0 / 5760 * G.TIMESTEP_AGENT; //4320 mim = 3 days
    public int aOB_DEATH;// = (int) (129600/G.TIMESTEP_AGENT); //72000 min = 12000 ts = 50 days
    public int unit_RESORPTION;// = (int) (1440 / G.TIMESTEP_AGENT); //1440 min = 240 ts = 1 day
    public int unit_FORMATION;//static = (int) (4320 / G.TIMESTEP_AGENT); //4320 min = 720 ts = 3 days

//    public double pOC_DEATH = 1.0 / 20160 * (MinToHour* TIMESTEP_AGENT); //1/(5760 min) = 1/(960 ts) = 1/(4 days)
//    public int max_aOB_DEATH = (int) (216000/(MinToHour* TIMESTEP_AGENT)); //129600 min = 21600 ts = 90 days; 216000 min = 150 days
//    public double MAX_pOB_DIFF_RATE = 1.0 / 20160 * (MinToHour* TIMESTEP_AGENT); //1/(20160 min) = 1/(3360 ts) = 1/(14 days)
    public int pOB_DIFF = (int) (20160 / (MinToHour* TIMESTEP_AGENT)); //20160 min = 3360 ts = 14 days
    public int max_RANKL_on = (int) (20160/(MinToHour* TIMESTEP_AGENT)); //20160 min = 3360 ts = 14 days (same as pOC_DEATH)
    double MSC_radius = 40.0 / SPACESTEP; // 40 microns (4 cells) Changed from 20 to 10 because 20 was too far in some cases
    double protect_radius = 20.0 / SPACESTEP; // 40 microns (4 cells)
//    double OB_radius = 200.0 / SPACESTEP; //arbitrarily set radius to 200
    //double MM_radius = 40.0 / G.SPACESTEP;
    double TGFBthresh = (1.05)* TGFB_basalRate/(Math.abs(TGFB_decayRate)* maxTGFB);//was 0.01; now 5% increase from basal

    //PROBABILISTIC RATES (UNIT TIME = 1 HOUR)
    public double MAX_FUSION_RATE = 1.0 / 4320 * (MinToHour); //*G.TIMESTEP_AGENT); //1/(4320 min) = 1/(720 ts) = 1/(3 days)
    public double pOB_DEATH = 1.0 / 4320 * (MinToHour); //*G.TIMESTEP_AGENT); // 4320 = 3 days;1/(20160 min) = 1/(3360 ts) = 1/(14 days)


    //ARRAY LISTS
    ArrayList<BoneCell_2022May17> tempPOClist = new ArrayList<>();
    ArrayList<BoneCell_2022May17> tempAOClist = new ArrayList<>();
    ArrayList<BoneCell_2022May17> LINING_recursionList = new ArrayList<>();
    ArrayList<BoneCell_2022May17> pOCrecursionList = new ArrayList<>();
    ArrayList<BoneCell_2022May17> aOCrecursionList = new ArrayList<>();
//    ArrayList<BoneCell_2022May17> aOBrecursionList = new ArrayList<>();
    ArrayList<BoneCell_2022May17> FusionList = new ArrayList<>();
    ArrayList<BoneCell_2022May17> ResorbingAOC = new ArrayList<>();
//    ArrayList<Integer> ResorbList = new ArrayList<>();
//    ArrayList<Integer> OC_IndexList = new ArrayList<>();
//    ArrayList<Integer> Current_OC_IndexList = new ArrayList<>();
//    ArrayList<Integer> Future_OC_IndexList = new ArrayList<>();


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CELL METHODS//////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //1. CheckPOC:   check to see if there are 5 adjacent pOC that are touching bone so that fusion occurs//
    //2. seekRANKL:                                                                                       //
    //3. ResorbStep: determine which bone to resorb provided aOC remains connected to neighbor aOC        //
    //4. CellStep:                                                                                        //
    ////////////////////////////////////////////////////////////////////////////////////////////////////////


    public ArrayList<BoneCell_2022May17> CheckPOC(ArrayList<BoneCell_2022May17> pOCrecursionList) {
        do{
            int rSize = pOCrecursionList.size();
            for (int i = 0; i < rSize; i++) { //Number of cells to check von Neumann neighborhood
                if(FusionList.size()!=5 && !FusionList.contains(pOCrecursionList.get(i))){
                    FusionList.add(pOCrecursionList.get(i));
                    fOpt = pOCrecursionList.get(i).MapOccupiedHood(G.MHood);//MapOccupiedHood(G.VNHood2); //Number of cells in VN neighborhood
                    for (int j = 0; j < fOpt; j++) {
                        //if (G.GetAgent(G.VNHood2[j]).type == pOC && G.GetAgent(G.VNHood2[j]).BoneInHood() == true && !FusionList.contains(G.GetAgent(G.VNHood2[j]))) {
                        //if (G.GetAgent(G.MHood[j]).type == pOC && (G.GetAgent(G.MHood[j]).BoneInHood() == true || G.GetAgent(G.MHood[j]).LiningInHood() == true) && !FusionList.contains(G.GetAgent(G.MHood[j]))) {
                        if (G.GetAgent(G.MHood[j]).type == pOC && G.GetAgent(G.MHood[j]).RanklInHood() == true && !FusionList.contains(G.GetAgent(G.MHood[j]))) {
                            //BoneInHood() uses VNHood so not redundant
                            //tempPOClist.add(G.GetAgent(G.VNHood2[j]));
                            tempPOClist.add(G.GetAgent(G.MHood[j]));
                        }
                    }
                }
            } //after this I have a list of neighboring pOCs
            pOCrecursionList.clear();
            pOCrecursionList.addAll(tempPOClist);
            tempPOClist.clear();
            //CheckPOC(pOCrecursionList);
        } while(FusionList.size()!=5 && pOCrecursionList.size()!=0);
        return FusionList;
    }


    public ArrayList<BoneCell_2022May17> CheckAOC(ArrayList<BoneCell_2022May17> aOCrecursionList) {
        do{
            //ResorbingAOC is cell list
            int rSize = aOCrecursionList.size();
            for (int i = 0; i < rSize; i++) { //Number of cells to check Moore neighborhood
                if(ResorbingAOC.size()<5 && !ResorbingAOC.contains(aOCrecursionList.get(i))){
                    ResorbingAOC.add(aOCrecursionList.get(i));
                    fOpt = aOCrecursionList.get(i).MapOccupiedHood(G.MHood);//MapOccupiedHood(G.VNHood2); //Number of cells in VN neighborhood
                    for (int j = 0; j < fOpt; j++) {
                        //if (G.GetAgent(G.VNHood2[j]).type == aOC && G.GetAgent(G.VNHood2[j]).BoneInHood() == true && !FusionList.contains(G.GetAgent(G.VNHood2[j]))) {
                        if (G.GetAgent(G.MHood[j]).type == aOC && G.GetAgent(G.MHood[j]).aOC_ID == this.aOC_ID  && !ResorbingAOC.contains(G.GetAgent(G.MHood[j]))) {
                            //tempAOClist.add(G.GetAgent(G.VNHood2[j]));
                            tempAOClist.add(G.GetAgent(G.MHood[j]));
                        }
                    }
                }
            } //after this I have a list of neighboring aOCs
            aOCrecursionList.clear();
            aOCrecursionList.addAll(tempAOClist);
            tempAOClist.clear();
            //CheckAOC(aOCrecursionList);
        } while(ResorbingAOC.size()<5 && aOCrecursionList.size()!=0);
        return ResorbingAOC;
    }



    public ArrayList<BoneCell_2022May17> CollectLining(ArrayList<BoneCell_2022May17> LINING_recursionList) {
        int N=0;
        int N_aOB=0;
        int N_other=0;
        do{
            int rSize = LINING_recursionList.size();
            for (int i = 0; i < rSize; i++) { //Number of cells to check Moore neighborhood
                if(!G.LiningList.contains(LINING_recursionList.get(i))){
                    G.LiningList.add(LINING_recursionList.get(i));
                    int eOpt = LINING_recursionList.get(i).MapOccupiedHood(G.MHood); //Number of cells in Moore neighborhood
                    for (int j = 0; j < eOpt; j++) { //Determine if cells are Bone lining and are not already in LiningList
                        //MarrowInHood uses VNHood so not redundant
                        if (G.GetAgent(G.MHood[j]).type == BONE && G.GetAgent(G.MHood[j]).RANKL_on == false && G.GetAgent(G.MHood[j]).BuriedInBone()==false && G.GetAgent(G.MHood[j]).ErodedBoneInHood()==false && G.GetAgent(G.MHood[j]).LiningInHood()==false &&  (!G.LiningList.contains(G.GetAgent(G.MHood[j])))) { // && G.GetAgent(G.MHood[j]).MarrowInHood()
                            G.tempEventList.add(G.GetAgent(G.MHood[j]));
                        } else if (G.GetAgent(G.MHood[j]).type == BONE && G.GetAgent(G.MHood[j]).RANKL_on == false && G.GetAgent(G.MHood[j]).BuriedInBone()==false && G.GetAgent(G.MHood[j]).ErodedBoneInHood()==false && G.GetAgent(G.MHood[j]).LiningInHood()==true &&  (!G.LiningList.contains(G.GetAgent(G.MHood[j])))){ //&& G.GetAgent(G.MHood[j]).MarrowInHood()
                            //G.tempLiningList.add(G.GetAgent(G.MHood[j]).Isq());
                            G.LiningList.add(G.GetAgent(G.MHood[j]));
                            N++;
                        //} else if(G.GetAgent(G.MHood[j]).type == LINING && !G.tempLiningList.contains(G.GetAgent(G.MHood[j]).Isq())){
                            //G.tempLiningList.add(G.GetAgent(G.MHood[j]).Isq());
                            //N++;
                        } else if(G.GetAgent(G.MHood[j]).type == aOB && G.GetAgent(G.MHood[j]).BuriedInBone()==false && G.GetAgent(G.MHood[j])!=this && !G.tempOBList.contains(G.GetAgent(G.MHood[j]).Isq())){ //&& G.GetAgent(G.MHood[j]).MarrowInHood() //this=current aOB this is dying; MarrowInHood to determine if aOB still capable of bone formation
                            G.tempOBList.add(G.GetAgent(G.MHood[j]).Isq());
                            N_aOB++;
                        } else if(G.GetAgent(G.MHood[j]).type == aOC || G.GetAgent(G.MHood[j]).type == pOB){ // || G.GetAgent(G.MHood[j]).type == MSC){
                            N_other++;
                        }
                    }
                }
            }
            LINING_recursionList.clear();
            LINING_recursionList.addAll(G.tempEventList);
            G.tempEventList.clear();
            //G.tempLiningList.clear();
        } while(LINING_recursionList.size()!=0 && N_other==0);// && N<2);
        G.tempLiningList.clear();
        G.tempOBList.clear();
        if(N_other>0 || N_aOB>0){// || N<2){ //to test if this does better job of defining lining
            G.LiningList.clear();
        }
        return G.LiningList;
    }




    public int seekRANKL() {
        //Will need to be careful when RANKL is not defined
        //It's not defined on bone but that should not be an issue because a cell cannot move there
//        public int[] moveHood = VonNeumannHood(true); //Have option of no movement
        int neighbors = MapHood(G.moveHood); //includes self
        double RANKL_this = G.RANKL.Get(G.moveHood[0]); //or G.RANKL.Get(this.Xsq(), this.Ysq());
        double RANKL_right = G.RANKL.Get(G.moveHood[1]); //or G.RANKL.Get(this.Xsq() + 1, this.Ysq());
        double RANKL_left = G.RANKL.Get(G.moveHood[2]); //or G.RANKL.Get(this.Xsq() - 1, this.Ysq());
        double RANKL_up = G.RANKL.Get(G.moveHood[3]); //or G.RANKL.Get(this.Xsq(), this.Ysq() + 1);
        double RANKL_down = G.RANKL.Get(G.moveHood[4]);// or G.RANKL.Get(this.Xsq(), this.Ysq() - 1);
        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); //cumulative probabilities


        for (int i = 0; i < neighbors; i++) {
            if ((G.GetAgent(G.moveHood[i]) == null || G.GetAgent(G.moveHood[i])==G.GetAgent(this.Isq()) || G.GetAgent(G.moveHood[i]).type==MM || (G.GetAgent(G.moveHood[i]).type==pOB && G.exposedBone[G.ItoX(G.moveHood[i])][G.ItoY(G.moveHood[i])]==false)) && !aOC_aOB_InHood(G.moveHood[i])) {
                emptyHood.add(G.moveHood[i]); //add index to list

                switch (i) {
                    case 0: //no movement
                        double P0 = 1 - 4 * G.pOC_DiffCoef - G.pOC_TaxisCoef * G.maxRANKL * (RANKL_right + RANKL_left - 4 * RANKL_this + RANKL_up + RANKL_down);
                        if(P0<0){
                            P0=0;
                            //System.out.println("Warning: P0<0");
                        }
                        cProbArray.add(P0);
                        ProbSum+=P0;
                        break;
                    case 1: //right
                        double P2 = G.pOC_DiffCoef + (G.pOC_TaxisCoef * G.maxRANKL)/ 4 * (RANKL_right - RANKL_left);
                        if(P2<0){
                            P2=0;
                            //System.out.println("Warning: P2<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P2 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P2);
                        }
                        ProbSum+=P2;
                        break;
                    case 2: //left
                        double P1 = G.pOC_DiffCoef - (G.pOC_TaxisCoef * G.maxRANKL)/ 4 * (RANKL_right - RANKL_left);
                        if(P1<0){
                            P1=0;
                            //System.out.println("Warning: P1<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P1 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P1);
                        }
                        ProbSum+=P1;
                        break;
                    case 3: //up
                        double P4 = G.pOC_DiffCoef + (G.pOC_TaxisCoef * G.maxRANKL)/ 4 * (RANKL_up - RANKL_down);
                        if(P4<0){
                            P4=0;
                            //System.out.println("Warning: P4<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P4 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P4);
                        }
                        ProbSum+=P4;
                        break;
                    case 4: //down
                        double P3 = G.pOC_DiffCoef - (G.pOC_TaxisCoef * G.maxRANKL)/ 4 * (RANKL_up - RANKL_down);
                        if(P3<0){
                            P3=0;
                            //System.out.println("Warning: P3<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P3 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P3);
                        }
                        ProbSum+=P3;
                        break;
                }

            }
        }

        int i = 0;
        int moveToIndex=G.moveHood[0];//emptyHood.get(0); //If ProbSum=0, cells don't move.

        while (i < cProbArray.size()) {
            if (0 < rnum && rnum <= cProbArray.get(0)/ProbSum) {
                moveToIndex = emptyHood.get(0); //keep as 0 because only want to enter once
                break;
            } else if (cProbArray.get(i)/ProbSum < rnum && rnum <= cProbArray.get(i + 1)/ProbSum) {
                moveToIndex = emptyHood.get(i+1);
                break;
            }
            i++;
        }
        return moveToIndex;
    }

    public int seekTGFB(double DiffCoef, double TaxisCoef) {
        //Will need to be careful when TGFB is not defined
        //It's not defined on bone but that should not be an issue because a cell cannot move there
//        public int[] moveHood = VonNeumannHood(true); //Have option of no movement
        int neighbors = MapHood(G.moveHood); //includes self
        double TGFB_this = G.TGFB.Get(G.moveHood[0]); //or G.TGFB.Get(this.Xsq(), this.Ysq());
        double TGFB_right = G.TGFB.Get(G.moveHood[1]); //or G.TGFB.Get(this.Xsq() + 1, this.Ysq());
        double TGFB_left = G.TGFB.Get(G.moveHood[2]); //or G.TGFB.Get(this.Xsq() - 1, this.Ysq());
        double TGFB_up = G.TGFB.Get(G.moveHood[3]); //or G.TGFB.Get(this.Xsq(), this.Ysq() + 1);
        double TGFB_down = G.TGFB.Get(G.moveHood[4]);// or G.TGFB.Get(this.Xsq(), this.Ysq() - 1);
        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); //cumulative probabilities


        for (int i = 0; i < neighbors; i++) {
            if ((G.GetAgent(G.moveHood[i]) == null || G.GetAgent(G.moveHood[i])==G.GetAgent(this.Isq()) || (G.GetAgent(G.moveHood[i]).type==MM && this.type==pOB) || (G.GetAgent(G.moveHood[i]).type==MSC && this.type==pOB))  && !aOC_aOB_InHood(G.moveHood[i])) {
                emptyHood.add(G.moveHood[i]); //add index to list

                switch (i) {
                    case 0: //no movement
                        double P0 = 1 - 4 * DiffCoef - TaxisCoef * maxTGFB * (TGFB_right + TGFB_left - 4 * TGFB_this + TGFB_up + TGFB_down);
                        if(P0<0){
                            P0=0;
                            //System.out.println("Warning: P0<0");
                        }
                        cProbArray.add(P0);
                        ProbSum+=P0;
                        break;
                    case 1: //right
                        double P2 = DiffCoef + (TaxisCoef * maxTGFB)/ 4 * (TGFB_right - TGFB_left);
                        if(P2<0){
                            P2=0;
                            //System.out.println("Warning: P2<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P2 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P2);
                        }
                        ProbSum+=P2;
                        break;
                    case 2: //left
                        double P1 = DiffCoef - (TaxisCoef * maxTGFB)/ 4 * (TGFB_right - TGFB_left);
                        if(P1<0){
                            P1=0;
                            //System.out.println("Warning: P1<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P1 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P1);
                        }
                        ProbSum+=P1;
                        break;
                    case 3: //up
                        double P4 = DiffCoef + (TaxisCoef* maxTGFB)/ 4 * (TGFB_up - TGFB_down);
                        if(P4<0){
                            P4=0;
                            //System.out.println("Warning: P4<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P4 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P4);
                        }
                        ProbSum+=P4;
                        break;
                    case 4: //down
                        double P3 = DiffCoef - (TaxisCoef* maxTGFB)/ 4 * (TGFB_up - TGFB_down);
                        if(P3<0){
                            P3=0;
                            //System.out.println("Warning: P3<0");
                        }
                        if(cProbArray.size()>0) {
                            cProbArray.add(P3 + cProbArray.get(cProbArray.size() - 1));
                        } else {
                            cProbArray.add(P3);
                        }
                        ProbSum+=P3;
                        break;
                }

            }
        }

        int i = 0;
        int moveToIndex=G.moveHood[0]; //emptyHood.get(0); //If ProbSum=0, cells don't move.

        while (i < cProbArray.size()) {
            if (0 < rnum && rnum <= cProbArray.get(0)/ProbSum) {
                moveToIndex = emptyHood.get(0); //keep as 0 because only want to enter once
                break;
            } else if (cProbArray.get(i)/ProbSum < rnum && rnum <= cProbArray.get(i + 1)/ProbSum) {
                moveToIndex = emptyHood.get(i+1);
                break;
            }
            i++;
        }
        return moveToIndex;
    }


    public void CellTick(){
        if(type==aOC){
            this.aOCtick++;
            this.aOCage++;
            //System.out.println(aOCtick);
        } else if(type==aOB){// && G.exposedBone[this.Xsq()][this.Ysq()] == false){
            //this.aOBtick++;
            this.aOBage++;
            //if(G.exposedBone[this.Xsq()][this.Ysq()] == false){
            this.aOBtick++; //this keeps track of bone formation
            //}
        } else if(type==pOB && G.exposedBone[this.Xsq()][this.Ysq()]==true){
            this.pOBage++;
        } else if(type==LINING || type==BONE){
            this.liningAge++;
        }
    }

    void Init() {
        if (type == BONE || type==LINING) { //Include aOB because they will turn into bone
            G.count_BA++;
        }
    }

    void Kill() {
        if (type == BONE || type==LINING) {
            G.count_BA--;
            G.AllBoneList.remove(this);
        }
        this.Dispose();
    }


    public void MM_RecruitMSC (double MSC_radius, int time) {
        int[] MSC_recHood = MooreHood(false); //recruit MSC adjacent to MM
        int[] MM_Hood = CircleHood(false, MSC_radius); //make sure other MSC not within radius

        //int N_MM, N_MSC, N_empty;

        //      RECRUIT EXTRA MSC
        int filled_opt = MapOccupiedHood(MM_Hood);
        int rec_opt = MapEmptyHood(MSC_recHood);
        int N_MSC = 0;

        if(rec_opt>0) {
            for(int i=0; i<filled_opt; i++){
                if(G.GetAgent(MM_Hood[i]).type == MSC){
                    N_MSC++;
                }
            }
            if(N_MSC==0){
                int MSCi = G.rn.Int(rec_opt);
                G.NewAgentSQ(MSC_recHood[MSCi]).type = MSC;
                G.GetAgent(MSC_recHood[MSCi]).InitAttributes(time);
                G.MSC_List.add(G.GetAgent(MSC_recHood[MSCi]));
            }
        }
    }

    public void RecruitMSC (double MSC_radius, ArrayList<BoneCell_2022May17> FusionList, int time) {
        BoneCell_2022May17 MSCtoRemove;
        ArrayList<BoneCell_2022May17> Temp_MSC_List = new ArrayList<>();

        int MSCiToAdd;
        ArrayList<Integer> i_Empty_CircleHood = new ArrayList<>();
        int[] MSC_recHood = CircleHood(false, MSC_radius);
        int[] emptyHood = MooreHood(false);
        //Order aOC by index
        Collections.sort(FusionList, Comparator.comparingInt(h -> h.Isq()));
        //Get agent in middle of aOC and recruit 1 MSC within radius of aOC
        int N_Empty_CircleHood = FusionList.get(FusionList.size()/2).MapEmptyHood(MSC_recHood);
//

        //      RECRUIT MSC DURING EACH REMODELING EVENT
        if(N_Empty_CircleHood>0 && FusionList.get(FusionList.size()/2).CellInHood(CircleHood(false,MSC_radius),MSC)==false) {
            //First determine if there is some space around the empty spot
            for(int i=0;i<N_Empty_CircleHood;i++){
                if(G.MapEmptyHood(emptyHood,MSC_recHood[i])>0 && aOC_aOB_InHood(MSC_recHood[i])==false){
                    i_Empty_CircleHood.add(i);
                }
            }
            //If there is space around empty spot, randomly select index. Otherwise, randomly select index of empty spot.
            if(i_Empty_CircleHood.size()>0){ //index of MSC_recHood where empty spot surrounded by empty spots and no aOC/aOB in hood
                MSCiToAdd = i_Empty_CircleHood.get(G.rn.Int(i_Empty_CircleHood.size()));
            } else{
                MSCiToAdd = G.rn.Int(N_Empty_CircleHood); //index of MSC_recHood
            }
//            do {
//                MSCiToAdd = G.rn.Int(N_Empty_CircleHood);
//            } while (G.MapEmptyHood(emptyHood,MSC_recHood[MSCiToAdd])==0);// //MSC needs room to divide

            //MSC_List created at initiation of simulation. The idea is to add to it if need more MSC
            while (G.MSC_List.size() > 0) {
                //Randomly remove another MSC from grid to maintain constant population

                MSCtoRemove = G.MSC_List.get(G.rn.Int(G.MSC_List.size()));//.Dispose(); //TODO: ONLY DISPOSE IF MSC NOT IN 20 DIAMETER RADIUS!
                //if MSC near aOC or MM, cannot recruit
                if (MSCtoRemove.aOCInCircleHood() == true || MSCtoRemove.MMInHood() == true) {
                    G.MSC_List.remove(MSCtoRemove);
                    Temp_MSC_List.add(MSCtoRemove);
                //otherwise, remove one MSC
                } else if ((G.MSC_List.size()+Temp_MSC_List.size()) > G.BMSCpop){
                    G.MSC_List.remove(MSCtoRemove);
                    if(MYELOMA){
//                        MSCtoRemove.Dispose();
                        Temp_MSC_List.add(MSCtoRemove); //If don't want MSC removal
                    } else {
                        MSCtoRemove.Dispose();
                    }
                } else {
                    G.MSC_List.remove(MSCtoRemove);
                    if(MYELOMA){
//                        MSCtoRemove.Dispose();
                        Temp_MSC_List.add(MSCtoRemove); //If don't want MSC removal
                    } else {
                        MSCtoRemove.Dispose();
                        break;
                    }
                }
            }

            //MSCtoRemove.Dispose();
            //Add new agent after removing one to ensure that you don't remove the new MSC
            G.MSC_List.addAll(Temp_MSC_List);
            Temp_MSC_List.clear();
            G.NewAgentSQ(MSC_recHood[MSCiToAdd]).type = MSC;
            G.GetAgent(MSC_recHood[MSCiToAdd]).InitAttributes(time);
            //G.GetAgent(MSC_recHood[MSCiToAdd]).eventtime=time;
            G.MSC_List.add(G.GetAgent(MSC_recHood[MSCiToAdd]));
//            } else {
//            G.NewAgentSQ(MSC_recHood[MSCiToAdd]).type = MSC;
//            G.MSC_List.AddAgent(G.GetAgent(MSC_recHood[MSCiToAdd]));
        }


        //      RECRUIT MM DURING FIRST REMODELING EVENT
        if(MYELOMA && G.MMstart==false){
            int MMiToAdd;
            i_Empty_CircleHood.clear();
            N_Empty_CircleHood = FusionList.get(FusionList.size()/2).MapEmptyHood(MSC_recHood);
//
            if(N_Empty_CircleHood>0) {
                for(int i=0;i<N_Empty_CircleHood;i++){
                    if(G.MapEmptyHood(emptyHood,MSC_recHood[i])>0 && aOC_aOB_InHood(MSC_recHood[i])==false){ //if some space around empty spot
                        i_Empty_CircleHood.add(i);
                    }
                }

                if(i_Empty_CircleHood.size()>0){
                    MMiToAdd = i_Empty_CircleHood.get(G.rn.Int(i_Empty_CircleHood.size()));
                } else{
                    MMiToAdd = G.rn.Int(N_Empty_CircleHood); //index of MSC_recHood
                }
//                do {
//                    MMiToAdd = G.rn.Int(N_Empty_CircleHood);
//                } while (G.MapEmptyHood(emptyHood, MSC_recHood[MMiToAdd]) == 0);// //MM needs room to divide
                G.NewAgentSQ(MSC_recHood[MMiToAdd]).type = MM;
                G.GetAgent(MSC_recHood[MMiToAdd]).InitAttributes(time);
                G.MMstart=true;
            }
        }
//        G.MMstart=true;
    }

    public double Mineralization_Time(double TGFBval){
        //Mineralization Time will be given as fold change to basal steady state (TGFB_basalRate/(TGFB_decayRate*maxTGFB))
        double scalefactor = (19.2*Math.abs(TGFB_decayRate)*maxTGFB)/TGFB_basalRate;
        double Y0 = 5.055;
        double plateau = 0.5379;
        double k = 0.1136*scalefactor;
        double yval = (Y0-plateau)*Math.exp(-k*TGFBval)+plateau;
        double basaltime = ((3.0*4320.0) / (MinToHour* TIMESTEP_AGENT)); //4320 min = 720 ts = 3 days; //I think I multiplied by 3 so that average aOB lifespan = 90 days
        return basaltime/yval;
    }

    public double Prob_Divide (double TGFBval, double MaxDivRate,double halfmax) {
        //double Ts = G.TGFB_basalRate/(Math.abs(G.TGFB_decayRate)*G.maxTGFB);
        double n=2;
        //double halfmax=Math.sqrt(3.0)*TGFBthresh;//1.05*Ts;//0.3;//
        return MaxDivRate*(1/(1+Math.pow((halfmax/TGFBval),n)));
//        if (TGFBval >= 0 && TGFBval < 1) {
//            return MaxDivRate*(-1/Math.log10(0.1*TGFBval));
//            //return MaxDivRate*(-1/Math.log10(0.1*(TGFBval-G.TGFB_basalRate/(Math.abs(G.TGFB_decayRate)*G.maxTGFB))));
//        } else {
//            return MaxDivRate;
//        }
    }

    public double Prob_MM_Divide (double TGFBval, double MaxDivRate,double halfmax) {
        //double Ts = G.TGFB_basalRate/(Math.abs(G.TGFB_decayRate)*G.maxTGFB);
        //double halfmax=Ts;//1.05*Ts;//0.3;//
        double n=2;
        return MaxDivRate*(1/(1+Math.pow((halfmax/TGFBval),n)));
    }


    public double Prob_Fusion (double RANKLval) {
        double halfmax=0.01;
        return 1/(1+(halfmax/RANKLval));
    }


    public double Hill_Repressor (double half_max, double coefficient, double concentration) {
        return 1/(1+Math.pow((concentration/half_max),coefficient));
    }

    public double Linear (double yint, double x2, double y2, double concentration) {
        double m = (y2-yint)/(x2-0);
        return m*concentration+yint;
    }

    public int FormBone () {
        ArrayList<Integer> PreferredIndex = new ArrayList<>(); //index of free space connected to bone
//        ArrayList<Integer> PreferredIndex2 = new ArrayList<>(); //index of free space connected to bone
        ArrayList<Integer> AltIndex = new ArrayList<>(); //index of free space not connected to bone
        ArrayList<Integer> AnyIndex = new ArrayList<>(); //index of free space

        int options = MapEmptyHood(G.VNHood); //empty positions surrounding aOB
            for (int j = 0; j < options; j++) {
                //Determine if empty position is still connected to bone
                int opt1 = G.MapOccupiedHood(G.VNHood2, G.VNHood[j]); //bone surrounding empty position
                int nbone = 0;
                int npOB = 0;
                int naOB = 0;
                for (int k = 0; k < opt1; k++) {
                    if (G.GetAgent(G.VNHood2[k]).type == BONE || G.GetAgent(G.VNHood2[k]).type == LINING) {
                        nbone++;
                    } else if (G.GetAgent(G.VNHood2[k]).type == pOB) {// || G.GetAgent(G.VNHood2[k]).type == LINING) { //ADDED TO PREVENT BURIED LINING
                        npOB++;
                    } else if (G.GetAgent(G.VNHood2[k]).type == aOB && G.GetAgent(G.VNHood2[k]) != this) {
                        naOB++;
                    }
                }
                if (opt1 < 4 && nbone > 0 && npOB == 0 && naOB == 0){// && G.InitBoneList.contains(G.VNHood[j])) { //these conditions are so that pOB don't get buried //
                    PreferredIndex.add(G.VNHood[j]);
                } else if (opt1 < 4){// && G.InitBoneList.contains(G.VNHood[j])) {
                    AltIndex.add(G.VNHood[j]);
                } else {
                    AnyIndex.add(G.VNHood[j]);
                }

            }

        if (PreferredIndex.size() > 0) {
            return PreferredIndex.get(G.rn.Int(PreferredIndex.size()));
//        } else if (PreferredIndex2.size() > 0) {
//            return PreferredIndex2.get(G.rn.Int(PreferredIndex2.size()));
        } else if (AltIndex.size() > 0){

            ArrayList<Integer> vMarrow = new ArrayList<>(Collections.nCopies(AltIndex.size(), 0)); //keeps track of number empty spots around each altindex
            ArrayList<Integer> maxMarrow = new ArrayList<>(); //List of indices of vMarrow that correspond to maximum

            for (int k = 0; k < AltIndex.size(); k++){
                int opt2 = G.MapEmptyHood(G.VNHood,AltIndex.get(k));
                vMarrow.set(k,opt2);
            }

            for (int m = 0; m < vMarrow.size(); m++) {
                if (vMarrow.get(m) == Collections.max(vMarrow)) { //todo maybe add condition if max(vBone)>1
                    maxMarrow.add(m);
                }
            }

            //randomly pick one of the indices
//            return AltIndex.get(G.rn.Int(maxMarrow.size()));
            return AltIndex.get(maxMarrow.get(G.rn.Int(maxMarrow.size())));
        } else if (AnyIndex.size() > 0){
            return AnyIndex.get(G.rn.Int(AnyIndex.size()));
        } else {
            return this.Isq();
        }

    }

    public boolean BoneInHood () {
        int[] BHood = VonNeumannHood(false);
        int Nbone = 0;
        int options = MapOccupiedHood(BHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(BHood[i]).type == BONE){// || G.GetAgent(G.VNHood[i]).type == LINING) { //originally lining cells were not part of this statement
                //if (G.GetAgent(G.VNHood[i]).type == BONE) {
                Nbone++;
            }
        }

        return Nbone > 0;
    }


    public boolean BuriedInBone () {
        int[] BHood = VonNeumannHood(false);
        int Nbone = 0;
        int options = MapOccupiedHood(BHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(BHood[i]).type == BONE || G.GetAgent(BHood[i]).type == LINING || G.GetAgent(BHood[i]).type == aOB){
                //if (G.GetAgent(G.VNHood[i]).type == BONE) {
                Nbone++;
            }
        }

        return Nbone == 4;
    }

    public int aOBage_transfer() {
        int aOB_index = this.Isq();
        double val;
        double min = Double.MAX_VALUE;
        int list_size = G.aOB_List.size();
        if(list_size>1){
            for(int i=0; i<list_size; i++){
                val = Dist(G.aOB_List.get(i).Xsq(),G.aOB_List.get(i).Ysq());
                if(val<=min && G.aOB_List.get(i)!=this){
                    min = val;
                    aOB_index = G.aOB_List.get(i).Isq();
                    if(G.aOB_List.get(i).aOB_DEATH==0){ //aOB_DEATH=0 when aOC_depth=0, which could happen if aOC subunit adjacent to other aOC subunit so it doesn't resorb bone
                        return aOB_index; //exit this method if we get to this point. not sure why I included this case though...
                    }
                }

            }

        }
        return aOB_index; //returns index of closest aOB
    }

    public double NearestBone() { //Used in Write_MM_birth_death output file
        double val;
        double min = Double.MAX_VALUE;
        int list_size = G.AllBoneList.size();
        if(list_size>0){
            for(int i=0; i<list_size; i++){
                val = Dist(G.AllBoneList.get(i).Xsq(),G.AllBoneList.get(i).Ysq());
                if(val<=min){
                    min = val;
                }
            }
        }
        return min; //returns distance to the nearest bone
    }

    public boolean RanklInHood () {
        int[] RanklHood = MooreHood(false);//VonNeumannHood
        int Nrankl = 0;
        int options = MapOccupiedHood(RanklHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(RanklHood[i]).RANKL_on==true){// || G.GetAgent(G.VNHood[i]).type == LINING) { //originally lining cells were not part of this statement
                Nrankl++;
            }
        }

        return Nrankl > 0;
    }

    public boolean LiningToBone () {
        int[] LHood = MooreHood(false);
        int Nbone = 0;
        int options = MapOccupiedHood(LHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(LHood[i]).type == BONE || G.GetAgent(LHood[i]).type == LINING){// || G.GetAgent(G.VNHood[i]).type == LINING) { //originally lining cells were not part of this statement
                //if (G.GetAgent(G.VNHood[i]).type == BONE) {
                Nbone++;
            }
        }

        // && MarrowInHood()==false) {
        return Nbone == 8;
    }

    public boolean LiningInHood () {
        int[] LHood = MooreHood(false);//was VNHood
        int Nlining = 0;
        int options = MapOccupiedHood(LHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(LHood[i]).type == LINING){
                Nlining++;
            }
        }

        return Nlining > 0;
    }


    public boolean aOC_aOB_InHood (int index) {
        int[] OC_Hood = MooreHood(false);
        int NaOC = 0;
        int options = G.MapOccupiedHood(OC_Hood,index); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(OC_Hood[i]).type == aOC || G.GetAgent(OC_Hood[i]).type == aOB) {
                NaOC++;
            }
        }

        return NaOC > 0;
    }


    public boolean aOCInCircleHood () {
        int[] CHood = CircleHood(false, MSC_radius);
        int NaOC = 0;
        int options = MapOccupiedHood(CHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(CHood[i]).type == aOC) {
                NaOC++;
            }
        }

        return NaOC > 0;
    }


    public boolean MarrowInHood () {
        int[] MarrowHood = VonNeumannHood(false); //4 neighbors
        int options = MapEmptyHood(MarrowHood); //occupied positions surrounding pOB

        return options > 0;
    }

    public boolean MMInHood () {
        //int[] MMHood = VonNeumannHood(false); //4 neighbors
        int[] MMHood = CircleHood(false, G.MM_radius); //4 neighbors
        int Ncancer = 0;
        int options = MapOccupiedHood(MMHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(MMHood[i]).type == MM) {
                Ncancer++;
            }
        }

        return Ncancer > 0;
    }

    public boolean CellInHood (int[] hood, int cell_type) {
        //int[] MMHood = VonNeumannHood(false); //4 neighbors
        int[] cellHood = hood; //4 neighbors
        int Ncell = 0;
        int options = MapOccupiedHood(cellHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(cellHood[i]).type == cell_type) {
                Ncell++;
            }
        }

        return Ncell > 0;
    }

    public boolean MMInHood3 () {
        int[] MMHood = MooreHood(false); //4 neighbors
        int Ncancer = 0;
        int options = MapOccupiedHood(MMHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(MMHood[i]).type == MM) {
                Ncancer++;
            }
        }

        return Ncancer > 2;
    }

    public boolean pOBInCircleHood (double radius) {
        int[] CHood = CircleHood(false, radius);
        int NpOB = 0;
        int options = MapOccupiedHood(CHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            if (G.GetAgent(CHood[i]).type == MSC || G.GetAgent(CHood[i]).type == pOB) {
                NpOB++;
            }
        }

        return NpOB > 0;
    }

    public boolean ErodedBoneInHood () {
        int[] ErodedHood = VonNeumannHood(false); //4 neighbors
        int Nbone = 0;
        int options = MapEmptyHood(ErodedHood); //occupied positions surrounding pOB
        for (int i = 0; i < options; i++) {
            //if (G.InitBoneList.contains(ErodedHood[i])) {
            if (G.exposedBone[G.ItoX(ErodedHood[i])][G.ItoY(ErodedHood[i])]==true){
                Nbone++;
            }
        }

        return Nbone > 0;
    }

    public boolean Obstacles (ArrayList<BoneCell_2022May17> FusionList) {
        ArrayList<Integer> vBone = new ArrayList<>(Collections.nCopies(4, 0));
        ArrayList<Integer> vOther = new ArrayList<>(Collections.nCopies(4, 0));
        //int[] vOther = new int[4];
        for (int i = 0; i < FusionList.size(); i++) {
            int options = FusionList.get(i).MapHood(G.VNHood);
            for (int j = 0; j < options; j++) {
                if (G.GetAgent(G.VNHood[j]) != null && (G.GetAgent(G.VNHood[j]).type == BONE || G.GetAgent(G.VNHood[j]).type == LINING)) {
                    vBone.set(j, vBone.get(j) + 1);
                } else if (G.GetAgent(G.VNHood[j]) != null && G.GetAgent(G.VNHood[j]).type != pOC) {
                    //vOther[j]++;
                    vOther.set(j, vOther.get(j) + 1);
                }
            }
        }
        //all directions have obstacles
        return Collections.min(vOther) > 0;
    }

    public void CellStep(FileIO Write_pOBborn, FileIO Write_pOBdiff, FileIO Write_aOBdeath, FileIO Write_MM_birth_death, FileIO Write_R_MM_clone, int time, List<Integer> RMeventTimes, double [] Cell_Counts){

//        //Production of RANKL by LINING cells. TODO: only select cells where event occurs
//        if(G.initEventList.contains(G.GetAgent(Isq()))){
//            G.GetAgent(Isq()).type=BONE;
//            G.RANKL.Add(Isq(), G.productionRate/maxRANKL);
//        }

        ///////
        //pOC//
        ///////

        if(type==pOC && time!=eventtime) {

            pOCrecursionList.clear();
            FusionList.clear();

            if (BoneInHood() == true) {
                pOCrecursionList.add(this);
                FusionList = CheckPOC(pOCrecursionList);
            }

            //Check if pOC dies
            if (aOC_aOB_InHood(this.Isq())) { //this is to keep pOC from interfering with aOC or aOB
                //
                ////////////
                //pOC dies//
                ////////////
                //
                Dispose();
                //pOC is randomly placed to replace dead pOC and keep pOC count constant
                do {
                    pOCx = G.rn.Int(G.xDim);
                    pOCy = G.rn.Int(G.yDim);
                } while ((G.PopAt(pOCx, pOCy) > 0 || aOC_aOB_InHood(G.I(pOCx, pOCy))) && G.Pop() < (G.xDim * G.yDim));

                if(G.Pop() < (G.xDim * G.yDim)) { //eventually the grid fills up and no more space to replace pOC
                    G.NewAgentSQ(pOCx, pOCy).type = this.type;
                    G.GetAgent(pOCx, pOCy).InitAttributes(time);
                }
            } else if (FusionList.size() == 5 && Obstacles(FusionList)==false) { //prevent fusion when there are obstacles that would prevent resorption
                //
                //////////////
                //pOC fusion//
                //////////////
                //
                Collections.sort(FusionList, Comparator.comparingInt(h -> h.Isq()));
                //Get agent in middle of aOC
                int MidIndex = FusionList.get(FusionList.size() / 2).Isq();
                double pfusion;

                if (BORTEZOMIB && TREATMENT_ON && INDIRECT_EFFECT) {
                    pfusion = Hill_Repressor(0.05, 2, G.dose) * Prob_Fusion(G.RANKL.Get(MidIndex));
                } else {
                    pfusion = Prob_Fusion(G.RANKL.Get(MidIndex));
                }

                if (G.rn.Double() < ProbScale(MAX_FUSION_RATE,TIMESTEP_AGENT)){ //time-dependent probability of fusion

                    if (G.rn.Double() < pfusion){ //RANKL-dependent probability of fusion
                        G.aOC_ID_counter++;

                        if (BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
                            aOC_DEATH = (int) (Hill_Repressor(0.1, 1, G.dose) * G.boundedGaussian((20160 / (MinToHour* TIMESTEP_AGENT)), (10080 / (MinToHour* TIMESTEP_AGENT) / 3), (10080 / (MinToHour* TIMESTEP_AGENT)), (30240 / (MinToHour* TIMESTEP_AGENT))));
                        } else if (MYELOMA){
                            aOC_DEATH = (int) (G.aOC_scale*G.boundedGaussian((20160 / (MinToHour* TIMESTEP_AGENT)), (10080 / (MinToHour* TIMESTEP_AGENT) / 3), (10080 / (MinToHour* TIMESTEP_AGENT)), (30240 / (MinToHour* TIMESTEP_AGENT))));
                        } else {
                            aOC_DEATH = (int) G.boundedGaussian((20160 / (MinToHour* TIMESTEP_AGENT)), (10080 / (MinToHour* TIMESTEP_AGENT) / 3), (10080 / (MinToHour* TIMESTEP_AGENT)), (30240 / (MinToHour* TIMESTEP_AGENT)));
                        }

                        //Collect 5 lining cells that are expressing RANKL and turn off (RANKL no longer expressed because fusion occurred)
                        int options = G.GetAgent(MidIndex).MapOccupiedHood(G.MHood);
                        G.recursionList.clear();
                        G.initEventList.clear();
                        for (int j = 0; j < options; j++) {
                            if (G.GetAgent(G.MHood[j]).RANKL_on == true) {
                                G.recursionList.add(G.GetAgent(G.MHood[j]));
                            }
                        }
                        G.initEventList = G.EndRANKL(G.recursionList);
                        for (int i = 0; i < G.initEventList.size(); i++) {
                            G.GetAgent(G.initEventList.get(i).Isq()).RANKL_on = false;
                            G.GetAgent(G.initEventList.get(i).Isq()).MM_RANKL_on = false;
                            G.GetAgent(G.initEventList.get(i).Isq()).RANKLtimer = 0;
                            G.GetAgent(G.initEventList.get(i).Isq()).eventtime = time;
                        }
                        G.initEventList.clear();

                        //Convert list of 5 pOC to 1 aOC
                        for (int i = 0; i < FusionList.size(); i++) { //TODO: WHAT HAPPENS WHEN POC TURNS TO AOC? DO WE ITERATE EACH TWICE IN TIMESTEP?
                            FusionList.get(i).type = aOC;
                            FusionList.get(i).aOCtick = 0;
                            FusionList.get(i).aOCage = 0;
                            FusionList.get(i).eventtime = time; //this is to prevent other pOC from being re-defined during this timestep
                            FusionList.get(i).aOC_ID = G.aOC_ID_counter;
                            FusionList.get(i).TGFB_on = true;
                            FusionList.get(i).aOC_DEATH = aOC_DEATH;

                            //assign unit_resorption when aOC fuses
                            if (BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
                                FusionList.get(i).unit_RESORPTION = (int) (2.0 * (1440.0 / (MinToHour* TIMESTEP_AGENT))); //2880 min = 480 ts = 2 days
                            } else {
                                FusionList.get(i).unit_RESORPTION = (int) (1440.0 / (MinToHour* TIMESTEP_AGENT)); //1440 min = 240 ts = 1 day
                            }

                            //Replace pOC to maintain constant population
                            do {
                                pOCx = G.rn.Int(G.xDim);
                                pOCy = G.rn.Int(G.yDim);
                            } while ((G.PopAt(pOCx, pOCy) > 0 || aOC_aOB_InHood(G.I(pOCx, pOCy))) && G.Pop() < (G.xDim * G.yDim));

                            if(G.Pop() < (G.xDim * G.yDim)) {
                                G.NewAgentSQ(pOCx, pOCy).type = pOC; //Replace pOC
                                G.GetAgent(pOCx, pOCy).InitAttributes(time);
                            }
                        }
                        Write_pOBborn.Write("aOC" + "," + this.aOC_ID + "," + time + "," + G.RANKL.Get(MidIndex) + "," + G.TGFB.Get(MidIndex) + "\n");

                        RecruitMSC(MSC_radius, FusionList, time);
                        return;

                    } else { //turn RANKL off, no fusion occurs

                        //Collect RANKL_ON
                        int options = G.GetAgent(MidIndex).MapOccupiedHood(G.MHood);
                        G.recursionList.clear();
                        G.initEventList.clear();
                        for (int j = 0; j < options; j++) {
                            if (G.GetAgent(G.MHood[j]).RANKL_on == true) {
                                G.recursionList.add(G.GetAgent(G.MHood[j]));
                            }
                        }
                        G.initEventList = G.EndRANKL(G.recursionList);
                        for (int i = 0; i < G.initEventList.size(); i++) {
                            G.GetAgent(G.initEventList.get(i).Isq()).RANKL_on = false;
                            G.GetAgent(G.initEventList.get(i).Isq()).MM_RANKL_on = false;
                            G.GetAgent(G.initEventList.get(i).Isq()).RANKLtimer = 0;
                            G.GetAgent(G.initEventList.get(i).Isq()).eventtime = time;
                        }
                        G.initEventList.clear();
                    }
                }
            } else {
                //
                /////////////
                //pOC moves//
                /////////////
                //
                int pOC_moveToIndex = seekRANKL(); //Need to define explicitly since seekRANKL has random values
                if (G.GetAgent(pOC_moveToIndex)!=null && (G.GetAgent(pOC_moveToIndex).type==MM || G.GetAgent(pOC_moveToIndex).type==pOB)) {
                    SwapPosition(G.GetAgent(pOC_moveToIndex));
                } else {
                    MoveSQ(pOC_moveToIndex);
                }
            }
        }


        ///////
        //aOC//
        ///////

        if(type==aOC && time!=eventtime) {

            //aOC dies if it is older than aOC_DEATH (mean: 14 days)
            if (aOCage >= aOC_DEATH) {
                ////////////
                //aOC dies//
                ////////////
                int TotResorb = 0; //To sum total bone resorbed

                if (this.aOC_DONE==false) { //this is to record TotResorb from all aOC units in resorbedBone matrix, which is recorded in write_aOBdeath and used to be used in NeedBuild for aOBdeath, but no longer used. good for bone homeostasis plots in R
                    //Step 1: Collect aOC
                    aOCrecursionList.clear();
                    ResorbingAOC.clear();
                    aOCrecursionList.add(this);
                    ResorbingAOC = CheckAOC(aOCrecursionList); //ResorbingAOC is a cell list

                    for (int i = 0; i < ResorbingAOC.size(); i++) {
                        TotResorb += ResorbingAOC.get(i).resorbtick;
                    }
                    for (int j = 0; j < ResorbingAOC.size(); j++) {
                        ResorbingAOC.get(j).aOC_DONE = true;
                        G.resorbedBone[ResorbingAOC.get(j).Xsq()][ResorbingAOC.get(j).Ysq()] = TotResorb;
                    }
                }

                //Turn-on TGF-beta where aOC dies to keep pOB local to resorption pit
                    G.extraTGFB[this.Xsq()][this.Ysq()] = true;

                Write_aOBdeath.Write(time + "," + "aOC" + "," + this.aOCage + "," + this.aOC_ID + "," + this.resorbtick + "," + G.count_BA + "," + G.resorbedBone[this.Xsq()][this.Ysq()]  + "\n");
                G.exposedBone[this.Xsq()][this.Ysq()] = true;
                G.aOC_depth[this.Xsq()][this.Ysq()] = resorbtick;
                Dispose();
            } else if (aOCtick >= unit_RESORPTION){// && BoneInHood() == true) {
                ////////////////////
                //aOC resorbs bone//
                ////////////////////

                /////////////////////////////////////////////
                //GROUP-RESORPTION (NEW WAY OF DEFINING OC)//
                /////////////////////////////////////////////

                //Step 1: Collect aOC
                aOCrecursionList.clear();
                ResorbingAOC.clear();
                aOCrecursionList.add(this);
                ResorbingAOC = CheckAOC(aOCrecursionList);

                //Step 2: Determine direction of most bone
                //int[] vBone = new int[4];
                ArrayList<Integer> vBone = new ArrayList<>(Collections.nCopies(4, 0));
                ArrayList<Integer> vOther = new ArrayList<>(Collections.nCopies(4, 0));
                //int[] vOther = new int[4];
                for (int i = 0; i < ResorbingAOC.size(); i++) {
                    int options = ResorbingAOC.get(i).MapHood(G.VNHood);
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(G.VNHood[j]) != null && (G.GetAgent(G.VNHood[j]).type == BONE || G.GetAgent(G.VNHood[j]).type == LINING)) {
                            vBone.set(j, vBone.get(j) + 1);
                        } else if (G.GetAgent(G.VNHood[j]) != null && G.GetAgent(G.VNHood[j]).aOC_ID != ResorbingAOC.get(i).aOC_ID) {
                            //vOther[j]++;
                            vOther.set(j, vOther.get(j) + 1);
                        }
                    }
                }

                //Step 3: Remove direction that has obstacles
                for(int i=0; i<vBone.size(); i++){
                    if(vOther.get(i)>0){
                        vBone.set(i,0);
                    }
                }

                ArrayList<Integer> maxDir = new ArrayList<>(); //List of indices corresponding to direction of max bone
//                ArrayList<Integer> no_Obstacles = new ArrayList<>(); //List of indices corresponding to no obstacles
                int resorb_dir= 10;


                //Determine direction(s) of maximum bone
                if (Collections.max(vBone)>0) {
                    for (int i = 0; i < vBone.size(); i++) {
                        if (vBone.get(i) == Collections.max(vBone)) { //todo maybe add condition if max(vBone)>1
                            maxDir.add(i);
                        }
                    }
                    //Determine resorption direction
                    if (maxDir.size() == 1) { //only one max dir
                        resorb_dir = maxDir.get(0);
                    } else if (maxDir.size() > 1) {
                        //randomly pick one of the indices
                        resorb_dir = maxDir.get(G.rn.Int(maxDir.size()));
                    }
                } else { //turn of BDF if no bone to resorb
                        for(int k=0;k<ResorbingAOC.size();k++) {
                            ResorbingAOC.get(k).TGFB_on=false;
                       }
                }


                //Step 4: Resorb bone
                if (resorb_dir != 10) {
                    //for (int k = 0; k < ResorbingAOC.size(); k++) {
                    int k=0;
                    while(ResorbingAOC.size()>0){
                        //int options = ResorbingAOC.get(k).MapHood(G.VNHood);
                        ResorbingAOC.get(k).MapHood(G.VNHood); //This maps VNHood

                        if (G.GetAgent(G.VNHood[resorb_dir]) != null && (G.GetAgent(G.VNHood[resorb_dir]).type == BONE || G.GetAgent(G.VNHood[resorb_dir]).type == LINING)) {
                            //if (G.GetAgent(G.VNHood2[resorb_dir]) != null && (G.GetAgent(G.VNHood2[resorb_dir]).type == BONE || G.GetAgent(G.VNHood2[resorb_dir]).type == LINING)) {
                            G.GetAgent(G.VNHood[resorb_dir]).Kill();
                            ResorbingAOC.get(k).MoveSQ(G.VNHood[resorb_dir]);
                            G.GetAgent(G.VNHood[resorb_dir]).resorbtick = G.GetAgent(G.VNHood[resorb_dir]).resorbtick + 1;
                            G.GetAgent(G.VNHood[resorb_dir]).TGFB_on = true;
                            G.GetAgent(G.VNHood[resorb_dir]).eventtime = time;
                            G.GetAgent(G.VNHood[resorb_dir]).aOCtick = 0;

                            //assign unit_resorption after each unit resorbed
                            if (BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
                                G.GetAgent(G.VNHood[resorb_dir]).unit_RESORPTION = (int) (2.0*(1440.0 / (MinToHour* TIMESTEP_AGENT))); //2.0 //1440 min = 240 ts = 1 day
                            } else {
                                G.GetAgent(G.VNHood[resorb_dir]).unit_RESORPTION = (int) (1440.0 / (MinToHour* TIMESTEP_AGENT)); //1440 min = 240 ts = 1 day
                            }

                            ResorbingAOC.remove(k);
                            k = Math.min(Math.abs(k--), 0);

                        } else if (G.GetAgent(G.VNHood[resorb_dir])==null){
                            //Free space so movement but no resorption
                            ResorbingAOC.get(k).MoveSQ(G.VNHood[resorb_dir]);
                            //G.GetAgent(G.VNHood[resorb_dir]).resorbtick = G.GetAgent(G.VNHood[resorb_dir]).resorbtick+1;
                            G.GetAgent(G.VNHood[resorb_dir]).TGFB_on = false;
                            G.GetAgent(G.VNHood[resorb_dir]).eventtime = time;
                            G.GetAgent(G.VNHood[resorb_dir]).aOCtick = 0;

                            //assign unit_resorption after each unit resorbed
                            if (BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
                                G.GetAgent(G.VNHood[resorb_dir]).unit_RESORPTION = (int) (2.0*(1440.0 / (MinToHour* TIMESTEP_AGENT))); //2.0 //1440 min = 240 ts = 1 day
                            } else {
                                G.GetAgent(G.VNHood[resorb_dir]).unit_RESORPTION = (int) (1440.0 / (MinToHour* TIMESTEP_AGENT)); //1440 min = 240 ts = 1 day
                            }

                            ResorbingAOC.remove(k);
                            k=Math.min(Math.abs(k--),0);
                        } else {
                            k++; //if there is nowhere to move (for example, when 2 aOC subunits in row)
                        }
                    }
                }
            }
        }

        ///////
        //MSC//
        ///////
        if(type==MSC && time!=eventtime){

            //MSC proliferates if sufficient TGFB
            double pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_MSC_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh);

//                //////////////
//                //Remove MSC//
//                //////////////
//            if(MYELOMA==true && CellInHood(CircleHood(false,MSC_radius),MM)==false && aOCInCircleHood()==false && G.MSC_List.size() > G.BMSCpop){ //MMInHood()
//                G.MSC_List.remove(this);
//                Dispose();
//            } else if

            if (G.TGFB.Get(this.Isq()) >= TGFBthresh && G.rn.Double() < ProbScale(pdiv, TIMESTEP_AGENT) && MMInHood()==false){//G.rn.Double() < pdiv){ // ) {

                ///////////////
                //MSC divides//
                ///////////////
                int options=MapEmptyHood(G.MHood);
                if(options>0) {//>3){//todo should this be options>3 like pOB?
                    //System.out.println(G.TGFB.Get(this.Xsq(),this.Ysq()));
                    ArrayList<Integer> no_OC = new ArrayList<>();
                    for (int i = 0; i < options; i++) {
                        if (!aOC_aOB_InHood(G.MHood[i])) {
                            no_OC.add(i);
                        }
                    }
                    if (no_OC.size() > 0) {
//                    int rint1=G.rn.Int(options);
                        int rint1 = no_OC.get(G.rn.Int(no_OC.size()));
                        G.NewAgentSQ(G.MHood[rint1]).type = pOB;
                        G.GetAgent(G.MHood[rint1]).InitAttributes(time);
                        //G.GetAgent(G.MHood[rint1]).eventtime=time;
                        //G.GetAgent(G.MHood[rint1]).pOBtick=0;
                        G.pOB_ID_counter++;
                        G.GetAgent(G.MHood[rint1]).pOB_ID = G.pOB_ID_counter;
                        Write_pOBborn.Write("pOB" + "," + G.GetAgent(G.MHood[rint1]).pOB_ID + "," + time + "," + G.RANKL.Get(this.Isq()) + "," + G.TGFB.Get(this.Isq()) + "\n");
                    }
                }
            } else if (CellInHood(MooreHood(false),MM)==false) {
                /////////////
                //MSC moves//
                /////////////
                //If MSC doesn't proliferate, MSC moves
                MoveSQ(seekTGFB(G.MSC_DiffCoef,G.MSC_TaxisCoef)); //MSC_seekTGFB
//            int options = MapEmptyHood(G.moveHood);
//            if(options>0) {
//                MoveSQ(G.moveHood[G.rn.Int(options)]);
//            }
            }
        }

        ///////
        //pOB//
        ///////
        if(type==pOB && time!=eventtime){
            RANKL_on=false;//true;

            double rn_BirthDeath = G.rn.Double();

            //pOB proliferates if sufficient TGFB
            double pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh);
            double pdeath=pOB_DEATH;

            if ((G.exposedBone[this.Xsq()][this.Ysq()] == false  && rn_BirthDeath < ProbScale(pdeath,TIMESTEP_AGENT)) || pOBage >= 3*pOB_DIFF){ //rn_BirthDeath < pdeath //&& G.TGFB.Get(this.Isq()) < TGFBthresh) //(MarrowInHood() == false || G.TGFB.Get(this.Isq()) < TGFBthresh//aOCInHood()==false//BoneInHood()==false &&
                //
                ////////////
                //pOB dies//
                ////////////
                //
                if(pOBage>=2*pOB_DIFF) {//Not sure why I have this restriction here
                    G.extraTGFB[this.Xsq()][this.Ysq()] = false;
                    G.exposedBone[this.Xsq()][this.Ysq()] = false;
                    G.TGFBtimer[this.Xsq()][this.Ysq()] = 0;
                    G.aOC_depth[this.Xsq()][this.Ysq()] = 0;
                    G.resorbedBone[this.Xsq()][this.Ysq()] = 0;
                }
                this.RANKL_on=false;
                Dispose();
            } else if(G.exposedBone[this.Xsq()][this.Ysq()]==true && G.TGFB.Get(this.Isq())<TGFBthresh && MMInHood()==false){//BoneInHood()==true) {

                    //
                    //////////////////////
                    //pOB differentiates//
                    //////////////////////
                    //
                    //System.out.println(G.TGFB.Get(this.Isq()));
                    //
                    if(pOBtick < pOB_DIFF){
                        pOBtick++;
                    } else {
                        this.type = aOB;
                        this.NeedBuild = G.resorbedBone[this.Xsq()][this.Ysq()];
                        this.eventtime=time;
                        this.aOBage = 0;
                        this.aOBtick=0;
                        G.aOB_List.add(this);

                        //assign bone formation rate
                        if(BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
//                            this.unit_FORMATION = (int) Mineralization_Time(G.TGFB.Get(this.Isq()));
                            this.unit_FORMATION = (int) (Mineralization_Time(G.TGFB.Get(this.Isq()))/1.7);
                        } else if (BORTEZOMIB && TREATMENT_ON && INDIRECT_EFFECT){
                            this.unit_FORMATION = (int) (Mineralization_Time(G.TGFB.Get(this.Isq()))/Linear(1.0,1.0, 2.0, G.dose));
                        } else {
                            this.unit_FORMATION = (int) Mineralization_Time(G.TGFB.Get(this.Isq()));
                        }

                        this.aOB_DEATH = (int) ((3.0*4320.0) / (MinToHour* TIMESTEP_AGENT)) * (G.aOC_depth[this.Xsq()][this.Ysq()]); //4320 minutes = 3 days

                        G.aOC_depth[this.Xsq()][this.Ysq()] = 0;
                        G.resorbedBone[this.Xsq()][this.Ysq()] = 0;
                        G.TGFBtimer[this.Xsq()][this.Ysq()]=0;
                        G.extraTGFB[this.Xsq()][this.Ysq()]=false;
                        G.exposedBone[this.Xsq()][this.Ysq()]=false; //changed this to be done after aob decides to form bone.
                        Write_pOBdiff.Write(this.pOB_ID + "," + time + "," + this.pOBage + "," + G.TGFB.Get(this.Isq()) + "\n");
                        //return;
                    }
                } else if (G.TGFB.Get(this.Isq()) >= TGFBthresh && rn_BirthDeath < (ProbScale(pdeath,TIMESTEP_AGENT) + ProbScale(pdiv,TIMESTEP_AGENT))){// rn_BirthDeath < (pdeath + pdiv)){ //todo determine if I should keep tgfbthresh here
                    //
                    ///////////////
                    //pOB divides//
                    ///////////////
                    //
                    int options=MapEmptyHood(G.MHood);
                    //int[] qHood = MooreHood(false); //quiescent
                    int npOB=0;
                    //int nq=0; //number of quiescent


                    if(options>0) {//options>3// && npOB<9){//todo if this creates too many could require nPOB<3
                        //if(options>0 && npOB<3  && nq<5){
                        //System.out.println(G.TGFB.Get(this.Xsq(),this.Ysq()));
                        ArrayList<Integer> no_OC = new ArrayList<>();
                        for (int i = 0; i < options; i++) {
                            if (!aOC_aOB_InHood(G.MHood[i])) {
                                no_OC.add(i);
                            }
                        }
                        if (no_OC.size() > 0) {
//                        int rint2=G.rn.Int(options);
                            int rint2 = no_OC.get(G.rn.Int(no_OC.size()));
                            G.NewAgentSQ(G.MHood[rint2]).type = pOB;
                            G.GetAgent(G.MHood[rint2]).InitAttributes(time);
                            //G.GetAgent(G.MHood[rint2]).eventtime=time;
                            //G.GetAgent(G.MHood[rint2]).pOBtick=0;
                            G.pOB_ID_counter++;
                            G.GetAgent(G.MHood[rint2]).pOB_ID = G.pOB_ID_counter;
                            Write_pOBborn.Write("pOB" + "," + G.GetAgent(G.MHood[rint2]).pOB_ID + "," + time + "," + G.RANKL.Get(this.Isq()) + "," + G.TGFB.Get(this.Isq()) + "\n");
                        }
                    }
                } else if (G.exposedBone[this.Xsq()][this.Ysq()]==false || ErodedBoneInHood()==true){//(G.exposedBone[this.Xsq()][this.Ysq()]==true && G.extraTGFB[this.Xsq()][this.Ysq()]==true)){//if (BoneInHood()==false) {
                    //
                    /////////////
                    //pOB moves//
                    /////////////
                    //
//                    MoveSQ(seekTGFB(G.pOB_DiffCoef,G.pOB_TaxisCoef));
                    int pOB_moveToIndex = seekTGFB(G.pOB_DiffCoef,G.pOB_TaxisCoef); //Need to define explicitly since seekTGFB has random values
                    if (G.GetAgent(pOB_moveToIndex)!=null && (G.GetAgent(pOB_moveToIndex).type==MM || G.GetAgent(pOB_moveToIndex).type==MSC)) {
                        SwapPosition(G.GetAgent(pOB_moveToIndex));
                    } else {
                        MoveSQ(pOB_moveToIndex);
                    }
                }

//            if (G.exposedBone[this.Xsq()][this.Ysq()]==true && G.extraTGFB[this.Xsq()][this.Ysq()]==true){
//                G.extraTGFB[this.Xsq()][this.Ysq()]=false;
//            }
        }

        ///////
        //aOB//
        ///////
        if(type==aOB && time!=eventtime) {

            if (this.BuriedInBone()==true || aOBage >= aOB_DEATH) {

                //
                ////////////
                //aOB dies//
                ////////////
                //

                if(this.BuriedInBone()==true) {
                    if(aOBage<aOB_DEATH) {
                        if (G.GetAgent(aOBage_transfer())!=this) {
                            int AgeDiff = this.aOB_DEATH-this.aOBage;
                            G.GetAgent(aOBage_transfer()).aOB_DEATH +=  AgeDiff;
                        }
                    }
                    Write_aOBdeath.Write(time + "," + "aOB" + "," + this.aOBage + "," + this.pOB_ID + "," + this.formtick + "," + G.count_BA + "," + this.NeedBuild + "\n");
                    this.type=BONE;
                    this.Init();
                    this.eventtime=time;
                    G.aOB_List.remove(this);
                    G.AllBoneList.add(this);

                } else {
                    Write_aOBdeath.Write(time + "," + "aOB" + "," + this.aOBage + "," + this.pOB_ID + "," + this.formtick + "," + G.count_BA + "," + this.NeedBuild + "\n");
                    G.aOB_List.remove(this);
                    G.AllBoneList.remove(this);
                    Dispose();
                }
            } else if (aOBtick >= unit_FORMATION && MapEmptyHood(G.VNHood)>0 && G.exposedBone[this.Xsq()][this.Ysq()]==false){//BoneInHood()==true) {
                //
                //////////////////
                //aOB forms bone//
                //////////////////
                //

                // System.out.println(G.TGFB.Get(this.Xsq(),this.Ysq()));

                int aOB_index = this.Isq();
                //Move aOB
                MoveSQ(FormBone());
                if(this.Isq()!=aOB_index) { //bone doesn't form if aOB doesn't move
                    //Replace old position with bone
                    G.NewAgentSQ(aOB_index).type = BONE;
                    G.GetAgent(aOB_index).InitAttributes(time);
                    G.GetAgent(aOB_index).Init();
                    G.AllBoneList.add(G.GetAgent(aOB_index));
                    //G.GetAgent(aOB_index).eventtime=time;
                    this.aOBtick = 0; //Note to self: could use modulus with 1 timer instead of separate timers
                    this.formtick++;
                    if(BISPHOSPHONATE && TREATMENT_ON && G.dose!=0.0) {
//                        this.unit_FORMATION = (int) Mineralization_Time(G.TGFB.Get(this.Isq()));
                        this.unit_FORMATION = (int) (Mineralization_Time(G.TGFB.Get(this.Isq()))/1.7);
                    } else if (BORTEZOMIB && TREATMENT_ON && INDIRECT_EFFECT){
                        this.unit_FORMATION = (int) (Mineralization_Time(G.TGFB.Get(this.Isq()))/Linear(1.0,1.0, 2.0, G.dose));
                    } else{
                        this.unit_FORMATION = (int) Mineralization_Time(G.TGFB.Get(this.Isq()));
                    }
//                    if(G.exposedBone[this.Xsq()][this.Ysq()]==true){
//                        G.exposedBone[this.Xsq()][this.Ysq()]=false;
//                        G.extraTGFB[this.Xsq()][this.Ysq()]=false;
//                        G.TGFBtimer[this.Xsq()][this.Ysq()]=0;
//                        //this condition is in case aOB moves to position where there used to be aOC but no pOB filled spot so bone was still exposed.
//                        //if we don't set it to false, aOB cannot build bone. If pOB differentiates into aOB later at that spot, probably OK since they were supposed to before.
//                    }
                }
            } //else if (BoneInHood()==false) { //Not sure why I have aOB moving?!?
            //If aOB moves
            //MoveSQ(seekTGFB());
            //}
        }

        ////////
        //BONE//
        ////////

        if(type==BONE && RANKL_on==true && time!=eventtime) {
            if (RANKLtimer >= max_RANKL_on) {
                G.recursionList.clear();
                G.initEventList.clear();
                G.recursionList.add(this);
                G.initEventList = G.EndRANKL(G.recursionList);
                if (G.initEventList.size() == 5 && this.MM_RANKL_on==false) {
//                    int[] RMevents = new int[1];
//                    G.rn.RandomIS(RMevents, time, G.curI*G.TURNOVER_TIME);
                    int[] RMevents = new int[]{time+1};
                    List<Integer> newRMeventTimes = Arrays.stream(RMevents).boxed().collect(Collectors.toList()); //convert array to list
                    RMeventTimes.addAll(newRMeventTimes);
                }
                G.EndRemodelingEvent(G.initEventList,time);
            } else {
                RANKLtimer++;
            }
        } else if (type==BONE && RANKL_on==false && G.TGFB.Get(this.Isq())<TGFBthresh && time!=eventtime && BuriedInBone()==false && ErodedBoneInHood()==false){// MarrowInHood()==true){
                LINING_recursionList.clear();
                G.LiningList.clear();
                LINING_recursionList.add(this);
                G.LiningList = CollectLining(LINING_recursionList); //this is to define lining cells
                if(G.LiningList.size()>0){
                    for(int i=0;i<G.LiningList.size();i++){
                        G.LiningList.get(i).type=LINING;
                        G.LiningList.get(i).liningAge=0;
                        G.LiningList.get(i).eventtime=time;
                    }
                }
        }

        //could check if bone&touching marrow collectlining()

        //////////
        //LINING//
        //////////

        if(type==LINING && time!=eventtime) {
            if(LiningToBone()==true) {
                this.type = BONE;
//            }//single event
            } else if (RMeventTimes.contains(time) && this.liningAge>=G.TURNOVER_TIME/4 && BuriedInBone()==false) { //&& MarrowInHood()
                G.recursionList.clear();
                G.initEventList.clear();
                G.recursionList.add(this);
                G.initEventList = G.InitRANKL(G.recursionList);
                if (G.initEventList.size() == 5){// && G.SecondLayerBone(G.initEventList)==true) {
                    G.RemodelingEvent(G.initEventList, time);
                    //G.RemodelingEvent(this,time); //will need to recollect lining in lining list when all sites have been remodeled, or all consecutive lining <5.
                    G.last_event = time;
                    RMeventTimes.remove((Integer) time);
                }
//            }
            } else if(MYELOMA && MMInHood3() && BuriedInBone()==false){ //&& G.TGFB.Get(this.Isq())>TGFBthresh ){
                G.recursionList.clear();
                G.initEventList.clear();
                G.recursionList.add(this);
                G.initEventList = G.InitRANKL(G.recursionList);
                if (G.initEventList.size() == 5){// && G.SecondLayerBone(G.initEventList)==true) {
                    for(int i=0;i<G.initEventList.size();i++) {
                        G.GetAgent(G.initEventList.get(i).Isq()).MM_RANKL_on=true;
                    }
                        G.RemodelingEvent(G.initEventList,time);
                }
            }
        }

//        if (type==LINING && time!=eventtime && LiningToBone()==true){
//            this.type=BONE;
//        } else if(type==LINING && RMeventTimes.contains(time) && time!=eventtime){
//            //if(G.rn.Double() < G.TURNOVER_TIME){
//            //if(RANKL_on==false) {
//            G.RemodelingEvent(this,time); //will need to recollect lining in lining list when all sites have been remodeled, or all consecutive lining <5.
//            int eindex = RMeventTimes.indexOf(time);
//            RMeventTimes.remove(eindex);
//            //}
//            //}
//        }

        //System.out.println(Age());

        ///////////
        //MYELOMA//
        ///////////
        if(type==MM){
            color = (!this.RESISTANT) ? MM : BLACK;

            double rn_BirthDeath = G.rn.Double();

            double pdiv;
            double pdeath=0;

//            double pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh);

            //RESISTANT NEAR MSC/POB
            if(this.RESISTANT && G.pOBadv==true && this.pOBInCircleHood(protect_radius)==true){
//                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),2.0*G.MAX_pOB_DIVISION_RATE,G.Ts);
//                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(5.0/3.0)*G.Ts); //2R
//                pdiv = Prob_MM_Divide(G.Ts,G.MAX_RESISTANT_DIVISION_RATE,Math.sqrt(5.0/4.0)*G.Ts); //1R
                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_RESISTANT_DIVISION_RATE,Math.sqrt(5.0/4.0)*G.Ts); //1R

                //RESISTANT NO MSC/POB
            } else if(this.RESISTANT) {
//                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),0.5*G.MAX_pOB_DIVISION_RATE,G.Ts);
//                pdiv = Prob_MM_Divide(G.Ts,G.MAX_RESISTANT_DIVISION_RATE,Math.sqrt(5.0/4.0)*G.Ts);
                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_RESISTANT_DIVISION_RATE,Math.sqrt(5.0/4.0)*G.Ts);

                //SENSITIVE NEAR MSC/POB
            }else if(G.pOBadv==true && this.pOBInCircleHood(protect_radius)==true){
//                pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh); //same as pOB
//                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,G.Ts);
//                pdiv = Prob_MM_Divide(G.Ts,G.MAX_MM_DIVISION_RATE,G.Ts);
                pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_MM_DIVISION_RATE,G.Ts);

                //SENSITIVE NO MSC/POB
            } else {
//               pdiv = Prob_MM_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,G.Ts);
//               pdiv = Prob_Divide(G.TGFB_basalRate/(Math.abs(G.TGFB_decayRate)*G.maxTGFB),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh); //basalRate
//               pdiv = MAX_pOB_DIVISION_RATE; //maxRate
//               pdiv = Prob_Divide(G.Ts,G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh); //same as pOB
                pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh); //same as pOB

            }

            ///////////////////
            //MSC recruitment//
            ///////////////////

//          //int BONE = 0, pOB = 1, aOB = 2, pOC = 3, aOC = 4, MSC = 5, LINING = 6, MM = 7;
            if(G.MSC_List.size()< Math.round((Cell_Counts[7]+Cell_Counts[8])/50.0)){
                MM_RecruitMSC(MSC_radius, time);
            }

            //pdiv = Prob_Divide(G.TGFB.Get(this.Isq()),G.MAX_pOB_DIVISION_RATE,Math.sqrt(3.0)*TGFBthresh); //same as pOB

//            if(time%12000==0){
//                Write_MM_birth_death.Write("all" + "," +  time + "," + NearestBone() + "," + G.RANKL.Get(this.Isq()) + "," + G.TGFB.Get(this.Isq()) + "\n");
//            }

            //if ( ((TREATMENT_ON==false || G.dose==0.0 || this.RESISTANT) && rn_BirthDeath < G.MM_DEATH && G.BDFadv && G.TGFB.Get(this.Isq()) < TGFBthresh) || (!this.RESISTANT && BORTEZOMIB && TREATMENT_ON && G.dose!=0.0 && EMDR==false && rn_BirthDeath <(G.MM_DEATH+G.dose*G.MM_DEATH)) || (!this.RESISTANT && BORTEZOMIB && TREATMENT_ON && G.dose!=0.0 && EMDR==true && rn_BirthDeath <(G.MM_DEATH+G.dose*G.MM_DEATH) && G.pOBadv && pOBInCircleHood(protect_radius)==false && G.BDFadv && G.TGFB.Get(this.Isq()) < TGFBthresh) ){ //&& G.pOBadv && pOBInCircleHood(protect_radius)==false && G.BDFadv && G.TGFB.Get(this.Isq()) < TGFBthresh)){// //19.0 //EMDR && G.pOBadv && pOBInCircleHood(protect_radius)==false && G.BDFadv && G.TGFB.Get(this.Isq()) < TGFBthresh)){
            //(G.rn.Double() < G.MM_DEATH && G.BDFadv==false)


            ///////////
            //MM dies//
            ///////////

            if((!TREATMENT_ON || G.dose==0.0 || this.RESISTANT) && G.BDFadv && G.TGFB.Get(this.Isq()) < TGFBthresh) { // && rn_BirthDeath < G.MM_DEATH) {
                pdeath = G.MM_DEATH; //No BDF survival advantage (no tx)

            } else if((!TREATMENT_ON || G.dose==0.0 || this.RESISTANT) && G.BDFadv && G.TGFB.Get(this.Isq()) >= TGFBthresh){ // && rn_BirthDeath < G.MM_DEATH_BDF){
//                pdeath = G.MM_DEATH; //No BDF survival advantage (no tx)
                pdeath =  G.MM_DEATH_BDF; //BDF survival advantage (no tx)

            } else if (TREATMENT_ON && G.dose!=0.0 && !this.RESISTANT && BORTEZOMIB && !EMDR){ // && rn_BirthDeath <(G.MM_DEATH+G.dose*1.5*G.MM_DEATH)) { //G.MM_DEATH+G.dose*9.0*G.MM_DEATH
                pdeath = G.MM_DEATH+G.dose*G.MM_DEATH_BTZ_FACTOR*G.MM_DEATH; //No EMDR everywhere

            } else if (TREATMENT_ON && G.dose!=0.0 && !this.RESISTANT && BORTEZOMIB && EMDR && G.BDFadv && G.pOBadv && G.TGFB.Get(this.Isq()) < TGFBthresh && pOBInCircleHood(protect_radius)==false){ ////// && rn_BirthDeath <(G.MM_DEATH+G.dose*1.5*G.MM_DEATH)) { //G.MM_DEATH+G.dose*9.0*G.MM_DEATH
                pdeath = G.MM_DEATH+G.dose*G.MM_DEATH_BTZ_FACTOR*G.MM_DEATH; //Non-protected S Cells

            } else if (TREATMENT_ON && G.dose!=0.0 && !this.RESISTANT && BORTEZOMIB && EMDR && G.BDFadv && G.pOBadv && (G.TGFB.Get(this.Isq()) >= TGFBthresh || pOBInCircleHood(protect_radius))){ ////// &&  rn_BirthDeath < G.MM_DEATH) { //10.0/3.0*G.MM_DEATH
                pdeath = G.MM_EMDR_DEATH; //Protected S Cells
            }

            if(rn_BirthDeath < ProbScale(pdeath,TIMESTEP_AGENT)){//rn_BirthDeath < pdeath){
                //MM dies
                color = Util.WHITE;
                G.MM_Death_Indices.add(this.Isq());
//                Write_MM_birth_death.Write("death" + "," +  time + "," + NearestBone() + "," + G.RANKL.Get(this.Isq()) + "," + G.TGFB.Get(this.Isq()) + "\n");
                Dispose();

            } else if (rn_BirthDeath < (ProbScale(pdeath,TIMESTEP_AGENT) + ProbScale(pdiv,TIMESTEP_AGENT))){//rn_BirthDeath < (pdeath + pdiv)){

            //////////////
            //MM divides//
            //////////////

            int options=MapEmptyHood(G.MHood);

                if(options>0){//option>3){
                    ArrayList<Integer> no_OC = new ArrayList<>();
                    for (int i = 0; i < options; i++) {
                        if (!aOC_aOB_InHood(G.MHood[i])) { //This rule is to try to prevent cells from getting in way of resorption
                            no_OC.add(i);
                        }
                    }
                    if (no_OC.size() > 0) {
//                    int rint2=G.rn.Int(options);
                        int rint2 = no_OC.get(G.rn.Int(no_OC.size()));
                        G.NewAgentSQ(G.MHood[rint2]).type = MM;
                        G.GetAgent(G.MHood[rint2]).InitAttributes(time);
                        G.GetAgent(G.MHood[rint2]).RESISTANT = this.RESISTANT;
                        if (!this.RESISTANT  && G.rn.Double() < G.pmutate && G.Start_Time!=0){ //mutation
                            G.R_MM_ID_counter++;
                            G.GetAgent(G.MHood[rint2]).RESISTANT = true;
                            G.GetAgent(G.MHood[rint2]).R_MM_ID = G.R_MM_ID_counter;
                            Write_R_MM_clone.Write(time + "," +  TREATMENT_ON + "," + G.dose + "," + (G.TGFB.Get(this.Isq())>=TGFBthresh) + "," + pOBInCircleHood(protect_radius)  + "\n");
                        } else if (this.RESISTANT){
                            G.GetAgent(G.MHood[rint2]).R_MM_ID =this.R_MM_ID;
                        }

                        G.MM_Division_Indices.add(this.Isq());
//                        Write_MM_birth_death.Write("birth" + "," +  time + "," + NearestBone() + "," + G.RANKL.Get(this.Isq()) + "," + G.TGFB.Get(this.Isq()) + "\n");
                    }
                }

            } else if(CellInHood(MooreHood(false),MSC)==false) {

            /////////////
            //MM moves//
            /////////////

            MoveSQ(seekTGFB(0,G.MSC_TaxisCoef));//only towards resorption

            }

        }


    }

}
