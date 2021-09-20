/*******************************************************************************/
/* aloha_sim.c     Code for simulating pure ALOHA                              */
/* 02/03/20                                                                    */
/*                                                                             */
/*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define S_DEBUG_N
//#define DEFER_GEN
//#define DEFER_PWS_N

/*******************************************************************************/
/* User changeable constants                                                   */
/*                                                                             */
/*******************************************************************************/

#define     No_of_Nodes         (25)
#define     Cap_G               (3.75)
#define     small_g1            (Cap_G/No_of_Nodes)
#define     small_g2            (Cap_G/No_of_Nodes)
#define     small_g3            (Cap_G/No_of_Nodes)
#define     small_g4            (Cap_G/No_of_Nodes)
#define     small_g5            (Cap_G/No_of_Nodes)
#define     small_g6            (Cap_G/No_of_Nodes)
#define     small_g7            (Cap_G/No_of_Nodes)
#define     small_g8            (Cap_G/No_of_Nodes)
#define     small_g9            (Cap_G/No_of_Nodes)
#define     small_g10           (Cap_G/No_of_Nodes)
#define     small_g11           (Cap_G/No_of_Nodes)
#define     small_g12           (Cap_G/No_of_Nodes)
#define     small_g13           (Cap_G/No_of_Nodes)
#define     small_g14           (Cap_G/No_of_Nodes)
#define     small_g15           (Cap_G/No_of_Nodes)
#define     small_g16           (Cap_G/No_of_Nodes)
#define     small_g17            (Cap_G/No_of_Nodes)
#define     small_g18            (Cap_G/No_of_Nodes)
#define     small_g19            (Cap_G/No_of_Nodes)
#define     small_g20           (Cap_G/No_of_Nodes)
#define     small_g21           (Cap_G/No_of_Nodes)
#define     small_g22           (Cap_G/No_of_Nodes)
#define     small_g23           (Cap_G/No_of_Nodes)
#define     small_g24           (Cap_G/No_of_Nodes)
#define     small_g25           (Cap_G/No_of_Nodes)

#define     Capital_G           (0.6)           /* in Erlang */


//#define     P_WS_P_1            (1)           /* probability P_WS when P=1  */
#define     Defer_Duration      (100)            /* in multiple pf packets */

#define     packet_duration     (1)             /* in ms */
#define     Simulation_Duration (100)      /* in ms */
#define     iterations          (1)

#define     EVENT_POOL_SIZE     (46900000)

#define     TRUE                (1)
#define     FALSE               (0)

/*******************************************************************************/
/* structure definitions                                                       */
/*                                                                             */
/*******************************************************************************/

typedef struct _NodeInfo{         /* content of a node */

    int Node_ID;

    int generated;                /* number of packets generated */
    int dropped_deferred;         /* number of packets dropped or deferred for various reasons */
    int transmitted;              /* number of packets transmitted */
    int collided;                 /* number of packets collided with other nodes' Tx */
    int success;                  /* number of packets successfully Tx-mitted */

    int Tx_Status;                /* TRUE or FALSE */
    int Coll_Status;              /* TRUE or FALSE */

    struct _NodeInfo   *Next;
    struct _NodeInfo   *Prev;
} NodeInfo;

/* event types */

#define Ev_P_GEN   1
#define Ev_P_DEF   2
#define Ev_End_Tx  3

typedef struct _Event{                /* content of a simulation event */

    unsigned char       Event_type;   /* 1: packet generate; 2: end of a Tx */
    double              Event_time;   /* in ms */

    struct _NodeInfo    *node ;       /* related node */

    int                 used;         /* TRUE or FALSE */

    struct _Event       *Next;
    struct _Event       *Prev;
} Event;

/*******************************************************************************/
/* arrays, lists, and related auxiliary variables; need initialization         */
/*                                                                             */
/*******************************************************************************/

NodeInfo        *Node_List[No_of_Nodes];         /* array containing all nodes */

Event           *Event_Pool[EVENT_POOL_SIZE];    /* event pool that provides events */
int             ev_access_ptr;                   /* variable used for optimizing the pool management */


Event           *Event_List = (Event *)NULL;     /* sorted linked list maintaining simulation events */
Event           *last_inserted_event;            /* variable used for optimizing the event lit search */


NodeInfo        *Channel;                        /* linked list of nodes that are in Tx state */

/*******************************************************************************/
/*  variables that need initialization                                  */
/*                                                                             */
/*******************************************************************************/

double          global_time;   /* in ms */
int             Sim_Stop;
float           tm;

/*******************************************************************************/
/* Misc. global variables .......                                              */
/*                                                                             */
/*******************************************************************************/

int             pf_total_gen_packets; /* across the entire network */
int             pf_total_tx_success;
int             sc_n1, sc_n2, sc_n3, sc_n4, sc_n5, gen1, gen2, gen3, gen4;
int             new_g;
int             def, ind_def;
int             coll;
int             collision_count;
int             n_id;
float           small_g;
float           self_coll;
int             topology[No_of_Nodes][No_of_Nodes];
float           s_g[No_of_Nodes];
int             gen_vec[No_of_Nodes];
int             scs_vec[No_of_Nodes];
float           s1[iterations], s2[iterations], s3[iterations],s4[iterations], s5[iterations], s6[iterations],s7[iterations], s8[iterations], s9[iterations], s10[iterations],s11[iterations], s12[iterations], s13[iterations], s14[iterations],s15[iterations], s16[iterations],s17[iterations], s18[iterations], s19[iterations], s20[iterations],s21[iterations], s22[iterations], s23[iterations], s24[iterations],s25[iterations], S[iterations];
/*******************************************************************************/
/* function prototypes  ..                                                     */
/*                                                                             */
/*******************************************************************************/

int             Create_Nodes(void);
double          randZeroToOne();
int             Create_Ev_Pool (void);
Event           *Get_Event(void);
void            Return_Event(Event *);

void            Insert_Event(unsigned char, double, NodeInfo *);
Event           *Remove_Front_Event(void);
void            Remove_Event(Event *);
void            Print_Event(Event *);
void            Print_Event_List(void);

void            Packet_Generate_Defer_Action(Event *);
void            End_Transmission_Action(Event *);

void            Initiate_Event(void);
void            Simulate(void);

float          get_exp(float);

/*******************************************************************************/
/* function definitions ..                                                     */
/*                                                                             */
/*******************************************************************************/

/* This one returns a random no. following exponential dist. */
float get_exp(float mean_value)
{
    return (-mean_value * log ((double)rand()/(double)RAND_MAX)) ;
//        return (-mean_value * log ((double)rand()*(-mean_value)/(double)RAND_MAX)) ;
}

/* Create all nodes and store them in an array */
int  Create_Nodes(void)
{
    char        *fn="Create_Nodes";
    int         cnt;
    NodeInfo    *node;

    for (cnt = 0; cnt < No_of_Nodes; cnt++){

        node = (NodeInfo  *)malloc(sizeof(NodeInfo));
        if (!node){
            printf ("%s: Failed to malloc NodeInfo .. terminating ..\n", fn);
            exit (0);
        }

        node->Node_ID = cnt;

        node->generated = 0;
        node->dropped_deferred= 0;
        node->transmitted = 0;
        node->collided = 0;
        node->success = 0;

        node->Tx_Status = FALSE;
        node->Coll_Status = FALSE;

        node->Next = (NodeInfo *)NULL;
        node->Prev = (NodeInfo *)NULL;

        Node_List[cnt] = node;
    }
    return (1);
}

/*  create a pool of event objects that can be used during the simulation */
int Create_Ev_Pool (void)
{
    char        *fn="Create_Ev_Pool";
    int         cnt;
    Event       *event;

    for (cnt = 0; cnt < EVENT_POOL_SIZE; cnt++){

        event = (Event *)malloc(sizeof(Event));
        if (!event){
            printf ("%s: Failed to malloc Event .. terminating ..\n", fn);
            exit (0);
        }

        event->Event_type = 0;
        event->Event_time = 0.0;

        event->node = (NodeInfo *)NULL;

        event->used = FALSE;

        event->Next = (Event *)NULL;
        event->Prev = (Event *)NULL;

        Event_Pool[cnt] = event;
    }
    ev_access_ptr = 0;
    return (1);
}

/* scans the event pool and gets an unused event that can be used by the simulation kernel */
Event *Get_Event(void)
{
    char    *fn="Get_Event";
    int     scan_start_ptr = ev_access_ptr;

    while (Event_Pool[ev_access_ptr]->used == TRUE){

        ev_access_ptr = (ev_access_ptr + 1) % EVENT_POOL_SIZE;     /* circular pool really */
        if (scan_start_ptr == ev_access_ptr){                      /* one full circle search yielded nothing */
            printf("%s: No more event memory available .. terminating ...\n", fn);
            exit (0);
        }
    }

    /* ev_access_ptr points to an unused event block; note that ev_access_ptr is an index to an array */

    Event_Pool[ev_access_ptr]->used = TRUE;

    return (Event_Pool[ev_access_ptr]);
}

/* Returns an event block to the circular buffer pool */
void Return_Event(Event *event)
{

    event->Event_type = 0;
    event->Event_time = 0.0;

    event->node = (NodeInfo *)NULL;

    event->used = FALSE;   // this is the main returning action ..

    event->Next = (Event *)NULL;
    event->Prev = (Event *)NULL;
}

/* routine for inserting an event in the simulation event list */
void Insert_Event(unsigned char ev_type, double time, NodeInfo *node)
{
    char            *fn="Insert_Event";
    Event           *event, *temp_event;
    int             Place_Found;

    event = Get_Event();

    event->Event_type = ev_type;
    event->Event_time = time;

    event->node = node;

    event->Next = (Event *)0;
    event->Prev = (Event *)0;

    if (time < global_time){
        printf ("%s: .... time < Global time!! Terminating ..", fn);
        Print_Event(event);
        exit (0);
    }

    if (!Event_List){ /* event list is empty */
        Event_List = event;
    }
    else{  /* The list isn't empty */
        if (time < Event_List->Event_time){ /* insert at the front */
            event->Next = Event_List;
            Event_List->Prev = event;
            Event_List = event;
        }
        else{  /* insert within the list */
            if ((!last_inserted_event) ||
                ((last_inserted_event) && (time <= last_inserted_event->Event_time))){
                temp_event = Event_List;
                /* insert the event in left of last_inserted_event */
            }
            else{
                temp_event = last_inserted_event;
                /* insert the event in right of last_inserted_event */
            }

            Place_Found = 0;
            while ((temp_event->Next) && (!Place_Found) ){
                if (time >= temp_event->Next->Event_time){
                    temp_event = temp_event->Next;
                }
                else{
                    Place_Found = 1;
                }
            }
            event->Prev = temp_event;
            event->Next = temp_event->Next;
            if (temp_event->Next){
                temp_event->Next->Prev = event;
            }
            temp_event->Next = event;
        }
    }
    last_inserted_event = event;
}

/*  remove the front event from the list for calling the handler/firing action  */
Event  *Remove_Front_Event(void)
{
    char         *fn="Remove_Front_Event";
    Event        *event;

    if (!Event_List){
        printf ("%s: Event_list is already empty .. terminating..\n", fn);
        exit (0);
    }

    if (Event_List == last_inserted_event){
        last_inserted_event = (Event *)NULL;
    }

    event = Event_List;
    if (event->Next){
        event->Next->Prev = (Event *)NULL;
    }
    Event_List = event->Next;
    return event;
}

/* remove an event from any arbitrary location of the event list */
void Remove_Event(Event *event)
{
    if(event == Event_List)
        Event_List = event->Next;
    else{
        event->Prev->Next = event->Next;
        if (event->Next)
            event->Next->Prev = event->Prev;
    }
    Return_Event(event); /* return event to the event pool for further use */
}

/* Routine for printing a single event. */
void Print_Event(Event *event)
{
    char    *fn="Print_Event";

    printf("\n%s: Details of the event follows ...\n", fn);
    printf("event_type = %d\n", event->Event_type);
    printf("event_time = %20.20f\n", event->Event_time);
    printf("event belongs to node %d\n\n", event->node->Node_ID);
}

/* Routine for printing the entire event list. */
void Print_Event_List(void)
{
    char        *fn="Print_Event_List";
    Event       *event;
    int         ev_count = 0;

    event = Event_List;
    while (event){
        Print_Event(event);
        ev_count += 1;
        event = event->Next;
    }
    printf ("Event_Count = %d\n", ev_count);
}

/*  handler/action routine for packet generation at a given node */
void   Packet_Generate_Defer_Action(Event *event)
{
    char        *fn="Packet_Generate_Defer_Action";
    NodeInfo    *node;
    int         Collision = FALSE;
//    float      small_g = (double)Capital_G/(double)(No_of_Nodes);
//    printf("%f",small_g);
    n_id = event->node->Node_ID;


    small_g = s_g[n_id];


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        gen_vec[n_id]++;
        new_g++;

        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        return;
    }

    /* schedule an end-of-Tx event */
    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);

    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */


        /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */



    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tm1 = Channel->Node_ID;
        //printf("%d\n", tm2);

        if(topology[n_id][tm1]==1){
            Collision = TRUE;
        }

    }
    Channel = event->node;

    event->node->Tx_Status = TRUE;

    event->node->transmitted++;

    if (Collision){
        /* set collision flag for involved nodes */
        node = Channel;
        while (node){
            node->Coll_Status = TRUE;
            node = node->Next;
        }
    }

    return;

}


/*  handler/action routine for end of packet transmission */
void    End_Transmission_Action(Event *event)
{
    if (!event->node->Coll_Status){
        event->node->success++;
        pf_total_tx_success++;
        int nd = event->node->Node_ID;

        scs_vec[nd]++;
    }else{
        event->node->collided++;
        coll++;
    }

    /* remove the node from the channel */

    if (Channel == event->node){  // the first node in the list
        if (Channel->Next)
            Channel->Next->Prev = (NodeInfo *)NULL;
        Channel = Channel->Next;
    }
    else{
        if (!event->node->Next){  // the last node ...
            event->node->Prev->Next = (NodeInfo *)NULL;
        }else{  // the node is somewhere inside the list
            event->node->Next->Prev = event->node->Prev;
            event->node->Prev->Next = event->node->Next;

        }
    }

    event->node->Tx_Status = FALSE;
    event->node->Coll_Status = FALSE;

    event->node->Next = (NodeInfo *)NULL;
    event->node->Prev = (NodeInfo *)NULL;
}

/*   Routine to put an initial event for each node in order to get things going  */
void Initiate_Events()
{
    char            *fn="Initiate_Events";
    int             ii;

    for(ii = 0; ii < No_of_Nodes; ii++)
        Insert_Event(Ev_P_GEN, 0.1* ((double)rand()/(double)RAND_MAX), Node_List[ii]);
}

/* simulation kernel */
void Simulate()
{
    char            *fn="Simulate";
    Event           *event;

    while (!Sim_Stop){

        event = Remove_Front_Event();
        global_time = event->Event_time;

        switch(event->Event_type){

            case Ev_P_GEN:

                Packet_Generate_Defer_Action (event);
                break;

            case Ev_P_DEF:

                Packet_Generate_Defer_Action (event);
                break;

            case Ev_End_Tx:

                End_Transmission_Action (event);
                break;

            default:

                printf ("%s: Undefined event_type = %d!! terminating ..\n", fn, event->Event_type);
                exit (0);
        }
#ifdef S_DEBUG
        Print_Event(event);
#endif

        Return_Event(event);

        if (global_time >= Simulation_Duration)
            Sim_Stop = TRUE;
//            printf("Global Time : %d\n",global_time);

#ifdef NEVER
        if (pf_total_gen_packets >= 1000000)
            Sim_Stop = TRUE;
#endif

    }

    float act_g = ((double)(new_g-def)*(double)packet_duration)/(double)Simulation_Duration;
    self_coll = (float)ind_def/(float)new_g;

    printf("\n \nSimulation Ends\n\n");
    printf("Total Success : %d\n",pf_total_tx_success);
    printf("Total Generated Packets : %d\n",new_g);
    printf("Total Deferred Packets : %d\n",def);
    printf("Total Collided Packets : %d\n",coll);
    printf("Self collisions : %f\n",self_coll);
    printf("Network Load : %f\n",Capital_G);
    printf("Individual loads : %f\n",small_g);
    printf("Actual Network Load : %f\n",act_g);
    printf("Node 1 Load : %f\n",small_g1);
    printf("Node 2 Load : %f\n",small_g2);
    printf("Global Time : %f\n",global_time);
    printf("Simulation Time: %d\n",Simulation_Duration);
    printf("Packet Duration : %d\n",packet_duration);
    printf("Resulting Capital_S = %10.10f\n", ((double)pf_total_tx_success/(double)global_time)*(double)packet_duration);
    printf("Resulting Capital_S_1 = %10.10f\n", ((double)pf_total_tx_success/(double)Simulation_Duration)*(double)packet_duration);
//    printf("Theoretical S: %f\n",(act_g)*(expf(-act_g)));
    printf("Success ratio : %f\n",(act_g)*(double)pf_total_tx_success/(double)(new_g-def));
    printf("Node1 Success:%d\tNode2 Success:%d\tNode3 Success:%d\tNode4 Success:%d\n",scs_vec[0],scs_vec[1], scs_vec[2], scs_vec[3]);
    //printf("Node1 packets:%d\tNode2 packets:%d\n",gen_vec[0],gen_vec[1]);
    //printf("s1:%f\ts2:%f\n",((double)scs_vec[0]/(double)global_time)*(double)packet_duration,((double)scs_vec[1]/(double)global_time)*(double)packet_duration);
    //printf("s3:%f\ts4:%f\n",((double)scs_vec[2]/(double)global_time)*(double)packet_duration,((double)scs_vec[3]/(double)global_time)*(double)packet_duration);


    for (int k = 0; k<No_of_Nodes; k++){
        printf("s%d:%f\n",k+1,((double)scs_vec[k]/(double)global_time)*(double)packet_duration);
    }

}

/* Main function that initializes things and starts simulation and then prints stats .. */

double randZeroToOne()
{
    return rand()/(double)(RAND_MAX);
}


int main()
{
    char    *fn="main";


    FILE     *fptr;


    fptr = fopen("program_pa25n_100ms_aloha_015.txt","a");


#ifdef S_DEBUG
    printf ("I am here ... \n");
#endif

    Create_Nodes();
    Create_Ev_Pool();

    Event_List = (Event *)NULL;
    last_inserted_event = (Event *)NULL;
    ev_access_ptr = 0;

    Channel = (NodeInfo *)NULL;

    global_time = 1.0e-16;
    Sim_Stop = FALSE;

    Initiate_Events();

    /* Initialize stat. variables */

    pf_total_gen_packets = 0;
    pf_total_tx_success = 0;
    gen1=0;
    gen2=0;
    gen3=0;
    gen4=0;
    sc_n1 = 0;
    sc_n2 = 0;
    sc_n3 = 0;
    sc_n4 = 0;
    new_g = 0;
    def = 0;
    ind_def = 0;
    coll = 0;





    s_g[0]=small_g1;
    s_g[1]=small_g2;
    s_g[2]=small_g3;
    s_g[3]=small_g4;
    s_g[4]=small_g5;
    s_g[5]=small_g6;
    s_g[6]=small_g7;
    s_g[7]=small_g8;
    s_g[8]=small_g9;
    s_g[9]=small_g10;
    s_g[10]=small_g11;
    s_g[11]=small_g12;
    s_g[12]=small_g13;
    s_g[13]=small_g14;
    s_g[14]=small_g15;
    s_g[15]=small_g16;
    s_g[16]=small_g17;
    s_g[17]=small_g18;
    s_g[18]=small_g19;
    s_g[19]=small_g20;
    s_g[20]=small_g21;
    s_g[21]=small_g22;
    s_g[22]=small_g23;
    s_g[23]=small_g24;
    s_g[24]=small_g25;
    //s_g[4]=small_g5;



    for(int k=0;k<No_of_Nodes;k++){
        gen_vec[k]=0;
    }

    for(int k=0;k<No_of_Nodes;k++){
        scs_vec[k]=0;
    }



    for(int i=0;i<No_of_Nodes;i++){
        for(int j=0;j<No_of_Nodes;j++){
            topology[i][j]=1;
            if(i==j){
                topology[i][j]=0;
            }
        }
    }




    /*topology[0][7]=0;
    topology[0][8]=0;
    topology[0][12]=0;
    topology[0][11]=0;
    topology[0][14]=0;
    topology[0][13]=0;
    topology[0][16]=0;
    topology[0][17]=0;
    topology[0][18]=0;
    topology[0][19]=0;
    topology[0][22]=0;
    topology[0][23]=0;
    topology[1][8]=0;
    topology[1][9]=0;
    topology[1][12]=0;
    topology[1][10]=0;
    topology[1][13]=0;
    topology[1][14]=0;
    topology[1][15]=0;
    topology[1][17]=0;
    topology[1][18]=0;
    topology[1][19]=0;
    topology[1][23]=0;
    topology[1][24]=0;
    topology[2][5]=0;
    topology[2][9]=0;
    topology[2][10]=0;
    topology[2][11]=0;
    topology[2][13]=0;
    topology[2][14]=0;
    topology[2][15]=0;
    topology[2][16]=0;
    topology[2][18]=0;
    topology[2][19]=0;
    topology[2][20]=0;
    topology[2][24]=0;
    topology[3][6]=0;
    topology[3][5]=0;
    topology[3][10]=0;
    topology[3][11]=0;
    topology[3][12]=0;
    topology[3][14]=0;
    topology[3][15]=0;
    topology[3][16]=0;
    topology[3][17]=0;
    topology[3][19]=0;
    topology[3][20]=0;
    topology[3][21]=0;
    topology[4][6]=0;
    topology[4][7]=0;
    topology[4][10]=0;
    topology[4][11]=0;
    topology[4][12]=0;
    topology[4][13]=0;
    topology[4][15]=0;
    topology[4][16]=0;
    topology[4][17]=0;
    topology[4][18]=0;
    topology[4][19]=0;
    topology[4][21]=0;
    topology[4][22]=0;
    topology[5][2]=0;
    topology[5][3]=0;
    topology[5][12]=0;
    topology[5][13]=0;
    topology[5][16]=0;
    topology[5][17]=0;
    topology[5][18]=0;
    topology[5][19]=0;
    topology[5][21]=0;
    topology[5][22]=0;
    topology[5][23]=0;
    topology[5][24]=0;
	 topology[6][3]=0;
 	 topology[6][4]=0;
 	 topology[6][13]=0;
 	 topology[6][14]=0;
 	 topology[6][15]=0;
 	 topology[6][17]=0;
 	 topology[6][18]=0;
 	 topology[6][19]=0;
 	 topology[6][20]=0;
 	 topology[6][22]=0;
 	 topology[6][23]=0;
 	 topology[6][24]=0;
	 topology[7][0]=0;
 	 topology[7][4]=0;
 	 topology[7][10]=0;
 	 topology[7][14]=0;
 	 topology[7][15]=0;
 	 topology[7][16]=0;
 	 topology[7][18]=0;
 	 topology[7][19]=0;
 	 topology[7][20]=0;
 	 topology[7][21]=0;
 	 topology[7][23]=0;
 	 topology[7][24]=0;
 	 topology[8][0]=0;
 	 topology[8][1]=0;
 	 topology[8][10]=0;
 	 topology[8][11]=0;
 	 topology[8][15]=0;
 	 topology[8][16]=0;
 	 topology[8][17]=0;
 	 topology[8][19]=0;
 	 topology[8][20]=0;
 	 topology[8][21]=0;
 	 topology[8][22]=0;
 	 topology[8][24]=0;
	 topology[9][1]=0;
 	 topology[9][2]=0;
 	 topology[9][11]=0;
 	 topology[9][12]=0;
 	 topology[9][15]=0;
 	 topology[9][16]=0;
 	 topology[9][17]=0;
 	 topology[9][18]=0;
 	 topology[9][20]=0;
 	 topology[9][21]=0;
 	 topology[9][22]=0;
 	 topology[9][23]=0;
	 topology[10][1]=0;
 	 topology[10][2]=0;
 	 topology[10][3]=0;
 	 topology[10][4]=0;
 	 topology[10][7]=0;
 	 topology[10][8]=0;
 	 topology[10][17]=0;
 	 topology[10][18]=0;
 	 topology[10][21]=0;
 	 topology[10][22]=0;
 	 topology[10][23]=0;
 	 topology[10][24]=0;
 	 topology[11][0]=0;
 	 topology[11][2]=0;
 	 topology[11][3]=0;
 	 topology[11][4]=0;
 	 topology[11][8]=0;
 	 topology[11][9]=0;
 	 topology[11][18]=0;
 	 topology[11][19]=0;
 	 topology[11][20]=0;
 	 topology[11][22]=0;
 	 topology[11][23]=0;
 	 topology[11][24]=0;
	 topology[12][0]=0;
 	 topology[12][1]=0;
 	 topology[12][3]=0;
 	 topology[12][4]=0;
 	 topology[12][5]=0;
 	 topology[12][9]=0;
 	 topology[12][15]=0;
 	 topology[12][19]=0;
 	 topology[12][20]=0;
 	 topology[12][21]=0;
 	 topology[12][23]=0;
 	 topology[12][24]=0;
	 topology[13][0]=0;
 	 topology[13][1]=0;
 	 topology[13][2]=0;
 	 topology[13][4]=0;
 	 topology[13][5]=0;
 	 topology[13][6]=0;
 	 topology[13][15]=0;
 	 topology[13][16]=0;
 	 topology[13][20]=0;
 	 topology[13][21]=0;
 	 topology[13][22]=0;
 	 topology[13][24]=0;
	 topology[14][0]=0;
 	 topology[14][1]=0;
 	 topology[14][2]=0;
 	 topology[14][3]=0;
 	 topology[14][6]=0;
 	 topology[14][7]=0;
 	 topology[14][16]=0;
 	 topology[14][17]=0;
 	 topology[14][20]=0;
 	 topology[14][21]=0;
 	 topology[14][22]=0;
 	 topology[14][23]=0;
	 topology[15][1]=0;
 	 topology[15][2]=0;
 	 topology[15][3]=0;
 	 topology[15][4]=0;
 	 topology[15][6]=0;
 	 topology[15][7]=0;
 	 topology[15][8]=0;
 	 topology[15][9]=0;
 	 topology[15][12]=0;
 	 topology[15][13]=0;
 	 topology[15][22]=0;
 	 topology[15][23]=0;
 	 topology[16][0]=0;
 	 topology[16][2]=0;
 	 topology[16][3]=0;
 	 topology[16][4]=0;
 	 topology[16][5]=0;
 	 topology[16][7]=0;
 	 topology[16][8]=0;
 	 topology[16][9]=0;
 	 topology[16][13]=0;
 	 topology[16][14]=0;
 	 topology[16][23]=0;
 	 topology[16][24]=0;
	 topology[17][0]=0;
 	 topology[17][1]=0;
 	 topology[17][3]=0;
 	 topology[17][4]=0;
 	 topology[17][5]=0;
 	 topology[17][6]=0;
 	 topology[17][8]=0;
 	 topology[17][9]=0;
 	 topology[17][10]=0;
 	 topology[17][14]=0;
 	 topology[17][20]=0;
 	 topology[17][24]=0;
	 topology[18][0]=0;
 	 topology[18][1]=0;
 	 topology[18][2]=0;
 	 topology[18][4]=0;
 	 topology[18][5]=0;
 	 topology[18][6]=0;
 	 topology[18][7]=0;
 	 topology[18][9]=0;
 	 topology[18][10]=0;
 	 topology[18][11]=0;
 	 topology[18][20]=0;
 	 topology[18][21]=0;
 	 topology[19][0]=0;
 	 topology[19][1]=0;
 	 topology[19][2]=0;
 	 topology[19][3]=0;
 	 topology[19][5]=0;
 	 topology[19][6]=0;
 	 topology[19][7]=0;
 	 topology[19][8]=0;
 	 topology[19][11]=0;
 	 topology[19][12]=0;
 	 topology[19][21]=0;
 	 topology[19][22]=0;
	 topology[20][2]=0;
 	 topology[20][3]=0;
 	 topology[20][6]=0;
 	 topology[20][7]=0;
 	 topology[20][8]=0;
 	 topology[20][9]=0;
 	 topology[20][11]=0;
 	 topology[20][12]=0;
 	 topology[20][13]=0;
 	 topology[20][14]=0;
 	 topology[20][17]=0;
 	 topology[20][18]=0;
	 topology[21][3]=0;
 	 topology[21][4]=0;
 	 topology[21][5]=0;
 	 topology[21][7]=0;
 	 topology[21][8]=0;
 	 topology[21][9]=0;
 	 topology[21][10]=0;
 	 topology[21][12]=0;
 	 topology[21][13]=0;
 	 topology[21][14]=0;
 	 topology[21][18]=0;
 	 topology[21][19]=0;
	 topology[22][0]=0;
 	 topology[22][4]=0;
 	 topology[22][5]=0;
 	 topology[22][6]=0;
 	 topology[22][8]=0;
 	 topology[22][9]=0;
 	 topology[22][10]=0;
 	 topology[22][11]=0;
 	 topology[22][13]=0;
 	 topology[22][14]=0;
 	 topology[22][15]=0;
 	 topology[22][19]=0;

	 topology[23][0]=0;
 	 topology[23][1]=0;
 	 topology[23][5]=0;
 	 topology[23][6]=0;
 	 topology[23][7]=0;
 	 topology[23][9]=0;
 	 topology[23][10]=0;
 	 topology[23][11]=0;
 	 topology[23][12]=0;
 	 topology[23][14]=0;
 	 topology[23][15]=0;
 	 topology[23][16]=0;
	 topology[24][1]=0;
 	 topology[24][2]=0;
 	 topology[24][5]=0;
 	 topology[24][6]=0;
 	 topology[24][7]=0;
 	 topology[24][8]=0;
 	 topology[24][10]=0;
 	 topology[24][11]=0;
 	 topology[24][12]=0;
 	 topology[24][13]=0;
 	 topology[24][16]=0;
 	 topology[24][17]=0;

*/


/*
    //topology[0][7]=0;
    topology[0][8]=0;
    //topology[0][12]=0;
    //topology[0][11]=0;
    topology[0][14]=0;
    topology[0][13]=0;
    topology[0][16]=0;
    topology[0][17]=0;
    //topology[0][18]=0;
    //topology[0][19]=0;
    topology[0][22]=0;
    //topology[0][23]=0;
    //topology[1][8]=0;
    topology[1][9]=0;
    //topology[1][12]=0;
    topology[1][10]=0;
    //topology[1][13]=0;
    topology[1][14]=0;
    //topology[1][15]=0;
    topology[1][17]=0;
    topology[1][18]=0;
    //topology[1][19]=0;
    topology[1][23]=0;
    //topology[1][24]=0;
    topology[2][5]=0;
    //topology[2][9]=0;
    topology[2][10]=0;
    topology[2][11]=0;
    //topology[2][13]=0;
    //topology[2][14]=0;
    //topology[2][15]=0;
    //topology[2][16]=0;
    topology[2][18]=0;
    topology[2][19]=0;
    //topology[2][20]=0;
    topology[2][24]=0;
    topology[3][6]=0;
    //topology[3][5]=0;
    //topology[3][10]=0;
    topology[3][11]=0;
    topology[3][12]=0;
    //topology[3][14]=0;
    topology[3][15]=0;
    //topology[3][16]=0;
    //topology[3][17]=0;
    topology[3][19]=0;
    topology[3][20]=0;
    //topology[3][21]=0;
    //topology[4][6]=0;
    topology[4][7]=0;
    //topology[4][10]=0;
    //topology[4][11]=0;
    topology[4][12]=0;
    topology[4][13]=0;
    topology[4][15]=0;
    topology[4][16]=0;
    //topology[4][17]=0;
    //topology[4][18]=0;
    topology[4][19]=0;
    topology[4][21]=0;
    //topology[4][22]=0;
    topology[5][2]=0;
    //topology[5][3]=0;
    //topology[5][12]=0;
    topology[5][13]=0;
    //topology[5][16]=0;
    //topology[5][17]=0;
    topology[5][18]=0;
    topology[5][19]=0;
    topology[5][21]=0;
    topology[5][22]=0;
    //topology[5][23]=0;
    //topology[5][24]=0;
	 topology[6][3]=0;
 	 //topology[6][4]=0;
 	 //topology[6][13]=0;
 	 topology[6][14]=0;
 	 topology[6][15]=0;
 	 //topology[6][17]=0;
 	 //topology[6][18]=0;
 	 topology[6][19]=0;
 	 //topology[6][20]=0;
 	 topology[6][22]=0;
 	 topology[6][23]=0;
 	 //topology[6][24]=0;
	 //topology[7][0]=0;
 	 topology[7][4]=0;
 	 topology[7][10]=0;
 	 //topology[7][14]=0;
 	 topology[7][15]=0;
 	 topology[7][16]=0;
 	 //topology[7][18]=0;
 	 //topology[7][19]=0;
 	 //topology[7][20]=0;
 	 //topology[7][21]=0;
 	 topology[7][23]=0;
 	 topology[7][24]=0;
 	 topology[8][0]=0;
 	 //topology[8][1]=0;
 	 //topology[8][10]=0;
 	 topology[8][11]=0;
 	 //topology[8][15]=0;
 	 topology[8][16]=0;
 	 topology[8][17]=0;
 	 //topology[8][19]=0;
 	 topology[8][20]=0;
 	 //topology[8][21]=0;
 	 //topology[8][22]=0;
 	 topology[8][24]=0;
	 topology[9][1]=0;
 	 //topology[9][2]=0;
 	 //topology[9][11]=0;
 	 topology[9][12]=0;
 	 //topology[9][15]=0;
 	 //topology[9][16]=0;
 	 topology[9][17]=0;
 	 topology[9][18]=0;
 	 topology[9][20]=0;
 	 topology[9][21]=0;
 	 //topology[9][22]=0;
 	 //topology[9][23]=0;
	 topology[10][1]=0;
 	 topology[10][2]=0;
 	 //topology[10][3]=0;
 	 //topology[10][4]=0;
 	 topology[10][7]=0;
 	 //topology[10][8]=0;
 	 //topology[10][17]=0;
 	 topology[10][18]=0;
 	 //topology[10][21]=0;
 	 //topology[10][22]=0;
 	 topology[10][23]=0;
 	 topology[10][24]=0;
 	 //topology[11][0]=0;
 	 topology[11][2]=0;
 	 topology[11][3]=0;
 	 //topology[11][4]=0;
 	 topology[11][8]=0;
 	 //topology[11][9]=0;
 	 //topology[11][18]=0;
 	 topology[11][19]=0;
 	 topology[11][20]=0;
 	 //topology[11][22]=0;
 	 //topology[11][23]=0;
 	 topology[11][24]=0;
	 //topology[12][0]=0;
 	 //topology[12][1]=0;
 	 topology[12][3]=0;
 	 topology[12][4]=0;
 	 //topology[12][5]=0;
 	 topology[12][9]=0;
 	 topology[12][15]=0;
 	 //topology[12][19]=0;
 	 topology[12][20]=0;
 	 topology[12][21]=0;
 	 //topology[12][23]=0;
 	 //topology[12][24]=0;
	 topology[13][0]=0;
 	 //topology[13][1]=0;
 	 //topology[13][2]=0;
 	 topology[13][4]=0;
 	 topology[13][5]=0;
 	 //topology[13][6]=0;
 	 //topology[13][15]=0;
 	 topology[13][16]=0;
 	 //topology[13][20]=0;
 	 topology[13][21]=0;
 	 topology[13][22]=0;
 	 //topology[13][24]=0;
	 topology[14][0]=0;
 	 topology[14][1]=0;
 	 //topology[14][2]=0;
 	 //topology[14][3]=0;
 	 topology[14][6]=0;
 	 //topology[14][7]=0;
 	 //topology[14][16]=0;
 	 topology[14][17]=0;
 	 //topology[14][20]=0;
 	 //topology[14][21]=0;
 	 topology[14][22]=0;
 	 topology[14][23]=0;
	 //topology[15][1]=0;
 	 //topology[15][2]=0;
 	 topology[15][3]=0;
 	 topology[15][4]=0;
 	 topology[15][6]=0;
 	 topology[15][7]=0;
 	 //topology[15][8]=0;
 	 //topology[15][9]=0;
 	 topology[15][12]=0;
 	 //topology[15][13]=0;
 	 //topology[15][22]=0;
 	 topology[15][23]=0;
 	 topology[16][0]=0;
 	 //topology[16][2]=0;
 	 //topology[16][3]=0;
 	 topology[16][4]=0;
 	 //topology[16][5]=0;
 	 topology[16][7]=0;
 	 topology[16][8]=0;
 	 //topology[16][9]=0;
 	 topology[16][13]=0;
 	 //topology[16][14]=0;
 	 //topology[16][23]=0;
 	 topology[16][24]=0;
	 topology[17][0]=0;
 	 topology[17][1]=0;
 	 //topology[17][3]=0;
 	 //topology[17][4]=0;
 	 //topology[17][5]=0;
 	 //topology[17][6]=0;
 	 topology[17][8]=0;
 	 topology[17][9]=0;
 	 //topology[17][10]=0;
 	 topology[17][14]=0;
 	 topology[17][20]=0;
 	 //topology[17][24]=0;
	 //topology[18][0]=0;
 	 topology[18][1]=0;
 	 topology[18][2]=0;
 	 //topology[18][4]=0;
 	 topology[18][5]=0;
 	 //topology[18][6]=0;
 	 //topology[18][7]=0;
 	 topology[18][9]=0;
 	 topology[18][10]=0;
 	 //topology[18][11]=0;
 	 //topology[18][20]=0;
 	 topology[18][21]=0;
 	 //topology[19][0]=0;
 	 //topology[19][1]=0;
 	 topology[19][2]=0;
 	 topology[19][3]=0;
 	 topology[19][5]=0;
 	 topology[19][6]=0;
 	 //topology[19][7]=0;
 	 //topology[19][8]=0;
 	 topology[19][11]=0;
 	 //topology[19][12]=0;
 	 //topology[19][21]=0;
 	 topology[19][22]=0;
	 //topology[20][2]=0;
 	 topology[20][3]=0;
 	 //topology[20][6]=0;
 	 //topology[20][7]=0;
 	 topology[20][8]=0;
 	 topology[20][9]=0;
 	 topology[20][11]=0;
 	 topology[20][12]=0;
 	 //topology[20][13]=0;
 	 //topology[20][14]=0;
 	 topology[20][17]=0;
 	 //topology[20][18]=0;
	 //topology[21][3]=0;
 	 topology[21][4]=0;
 	 topology[21][5]=0;
 	 //topology[21][7]=0;
 	 //topology[21][8]=0;
 	 topology[21][9]=0;
 	 //topology[21][10]=0;
 	 topology[21][12]=0;
 	 topology[21][13]=0;
 	 //topology[21][14]=0;
 	 topology[21][18]=0;
 	 //topology[21][19]=0;
	 topology[22][0]=0;
 	 //topology[22][4]=0;
 	 topology[22][5]=0;
 	 topology[22][6]=0;
 	 //topology[22][8]=0;
 	 //topology[22][9]=0;
 	 //topology[22][10]=0;
 	 //topology[22][11]=0;
 	 topology[22][13]=0;
 	 topology[22][14]=0;
 	 //topology[22][15]=0;
 	 topology[22][19]=0;

	 //topology[23][0]=0;
 	 topology[23][1]=0;
 	 //topology[23][5]=0;
 	 topology[23][6]=0;
 	 topology[23][7]=0;
 	 //topology[23][9]=0;
 	 topology[23][10]=0;
 	 //topology[23][11]=0;
 	 //topology[23][12]=0;
 	 topology[23][14]=0;
 	 topology[23][15]=0;
 	 //topology[23][16]=0;
	 //topology[24][1]=0;
 	 topology[24][2]=0;
 	 //topology[24][5]=0;
 	 //topology[24][6]=0;
 	 topology[24][7]=0;
 	 topology[24][8]=0;
 	 topology[24][10]=0;
 	 topology[24][11]=0;
 	 //topology[24][12]=0;
 	 //topology[24][13]=0;
 	 topology[24][16]=0;
 	 //topology[24][17]=0;
*/





    for(int i=0;i<iterations; i++){

    srand(time(0));
    Simulate();

    fprintf(fptr, "%f\n",((double)pf_total_tx_success/(double)global_time)*(double)packet_duration);

    s1[i]= ((double)scs_vec[0]/(double)global_time)*(double)packet_duration;
    s2[i]= ((double)scs_vec[1]/(double)global_time)*(double)packet_duration;
    s3[i]= ((double)scs_vec[2]/(double)global_time)*(double)packet_duration;
    s4[i]= ((double)scs_vec[3]/(double)global_time)*(double)packet_duration;
    s5[i]= ((double)scs_vec[4]/(double)global_time)*(double)packet_duration;
    s6[i]= ((double)scs_vec[5]/(double)global_time)*(double)packet_duration;
    s7[i]= ((double)scs_vec[6]/(double)global_time)*(double)packet_duration;
    s8[i]= ((double)scs_vec[7]/(double)global_time)*(double)packet_duration;
    s9[i]= ((double)scs_vec[8]/(double)global_time)*(double)packet_duration;
    s10[i]= ((double)scs_vec[9]/(double)global_time)*(double)packet_duration;
    s11[i]= ((double)scs_vec[10]/(double)global_time)*(double)packet_duration;
    s12[i]= ((double)scs_vec[11]/(double)global_time)*(double)packet_duration;
    s13[i]= ((double)scs_vec[12]/(double)global_time)*(double)packet_duration;
    s14[i]= ((double)scs_vec[13]/(double)global_time)*(double)packet_duration;
    s15[i]= ((double)scs_vec[14]/(double)global_time)*(double)packet_duration;
    s16[i]= ((double)scs_vec[15]/(double)global_time)*(double)packet_duration;
    s17[i]= ((double)scs_vec[16]/(double)global_time)*(double)packet_duration;
    s18[i]= ((double)scs_vec[17]/(double)global_time)*(double)packet_duration;
    s19[i]= ((double)scs_vec[18]/(double)global_time)*(double)packet_duration;
    s20[i]= ((double)scs_vec[19]/(double)global_time)*(double)packet_duration;
    s21[i]= ((double)scs_vec[20]/(double)global_time)*(double)packet_duration;
    s22[i]= ((double)scs_vec[21]/(double)global_time)*(double)packet_duration;
    s23[i]= ((double)scs_vec[22]/(double)global_time)*(double)packet_duration;
    s24[i]= ((double)scs_vec[23]/(double)global_time)*(double)packet_duration;
    s25[i]= ((double)scs_vec[24]/(double)global_time)*(double)packet_duration;
    S[i]= s1[i]+s2[i]+s3[i]+s4[i]+s5[i]+s6[i]+s7[i]+s8[i]+s9[i]+s10[i]+s11[i]+s12[i]+s13[i]+s14[i]+s15[i]+s16[i]+s17[i]+s18[i]+s19[i]+s20[i]+s21[i]+s22[i]+s23[i]+s24[i]+s25[i];
    }


    float avg_s1 =0.0;
    float avg_s2 =0.0;
    float avg_s3 =0.0;
    float avg_s4 =0.0;
    float avg_s5 =0.0;
    float avg_s6 =0.0;
    float avg_s7 =0.0;
    float avg_s8 =0.0;
    float avg_s9 =0.0;
    float avg_s10 =0.0;
    float avg_s11 =0.0;
    float avg_s12 =0.0;
    float avg_s13 =0.0;
    float avg_s14 =0.0;
    float avg_s15 =0.0;
    float avg_s16 =0.0;
    float avg_s17 =0.0;
    float avg_s18 =0.0;
    float avg_s19 =0.0;
    float avg_s20 =0.0;
    float avg_s21 =0.0;
    float avg_s22 =0.0;
    float avg_s23 =0.0;
    float avg_s24 =0.0;
    float avg_s25 =0.0;
    float avg_S = 0.0;

    for(int i=0;i<iterations; i++){
    avg_s1 = avg_s1+s1[i];
    avg_s2 = avg_s2+s2[i];
    avg_s3 = avg_s3+s3[i];
    avg_s4 = avg_s4+s4[i];
    avg_s5 = avg_s5+s5[i];
    avg_s6 = avg_s6+s6[i];
    avg_s7 = avg_s7+s7[i];
    avg_s8 = avg_s8+s8[i];
    avg_s9 = avg_s9+s9[i];
    avg_s10 = avg_s10+s10[i];
    avg_s11 = avg_s11+s11[i];
    avg_s12 = avg_s12+s12[i];
    avg_s13 = avg_s13+s13[i];
    avg_s14 = avg_s14+s14[i];
    avg_s15 = avg_s15+s15[i];
    avg_s16 = avg_s16+s16[i];
    avg_s17 = avg_s17+s17[i];
    avg_s18 = avg_s18+s18[i];
    avg_s19 = avg_s19+s19[i];
    avg_s20 = avg_s20+s20[i];
    avg_s21 = avg_s21+s21[i];
    avg_s22 = avg_s22+s22[i];
    avg_s23 = avg_s23+s23[i];
    avg_s24 = avg_s24+s24[i];
    avg_s25 = avg_s25+s25[i];
    avg_S =  avg_S+S[i];
    }



    printf("\n\n\n\n");
    printf("Avg. s1: %f\n", avg_s1/(float)iterations);
    printf("Avg. s2: %f\n", avg_s2/(float)iterations);
    printf("Avg. s3: %f\n", avg_s3/(float)iterations);
    printf("Avg. s4: %f\n", avg_s4/(float)iterations);
    printf("Avg. s5: %f\n", avg_s5/(float)iterations);
    printf("Avg. s6: %f\n", avg_s6/(float)iterations);
    printf("Avg. s7: %f\n", avg_s7/(float)iterations);
    printf("Avg. s8: %f\n", avg_s8/(float)iterations);
    printf("Avg. s9: %f\n", avg_s9/(float)iterations);
    printf("Avg. s10: %f\n", avg_s10/(float)iterations);
    printf("Avg. s11: %f\n", avg_s11/(float)iterations);
    printf("Avg. s12: %f\n", avg_s12/(float)iterations);
    printf("Avg. s13: %f\n", avg_s13/(float)iterations);
    printf("Avg. s14: %f\n", avg_s14/(float)iterations);
    printf("Avg. s15: %f\n", avg_s15/(float)iterations);
    printf("Avg. s16: %f\n", avg_s16/(float)iterations);
    printf("Avg. s17: %f\n", avg_s17/(float)iterations);
    printf("Avg. s18: %f\n", avg_s18/(float)iterations);
    printf("Avg. s19: %f\n", avg_s19/(float)iterations);
    printf("Avg. s20: %f\n", avg_s20/(float)iterations);
    printf("Avg. s21: %f\n", avg_s21/(float)iterations);
    printf("Avg. s22: %f\n", avg_s22/(float)iterations);
    printf("Avg. s23: %f\n", avg_s23/(float)iterations);
    printf("Avg. s24: %f\n", avg_s24/(float)iterations);
    printf("Avg. s25: %f\n", avg_s25/(float)iterations);
    printf("Avg. S: %f\n", avg_S/(float)iterations);
}
