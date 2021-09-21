/*******************************************************************************/
/* aloha_sim.c     Code for simulating pure ALOHA                              */
/* 02/16/20                                                                    */
/*                                                                             */
/*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define S_DEBUG_N
//#define DEFER_GEN
#define DEFER_PWS_N

/*******************************************************************************/
/* User changeable constants                                                   */
/*                                                                             */
/*******************************************************************************/

#define     No_of_Nodes         (12)
#define     small_g1            (0.3)
#define     small_g2            (0.3)
#define     small_g3            (0.3)
#define     small_g4            (0.3)
#define     small_g5            (0.3)
#define     small_g6            (0.3)
#define     small_g7            (0.3)
#define     small_g8            (0.3)
#define     small_g9            (0.3)
#define     small_g10            (0.3)
#define     small_g11            (0.3)
#define     small_g12            (0.3)
#define     packet_duration     (1)             /* in ms */
#define     Simulation_Duration (10000)        /* in ms */
#define     packets_per_epoch   (1500)
#define     EVENT_POOL_SIZE     (5500000)

#define     freq                (66)

#define     TRUE                (1)
#define     FALSE               (0)


#define     no_of_epochs        (30000)

#define     pkt_er_rate         (0.00)

#define     iterations          (8)

#define     alpha               (0.9)
#define     beta                (0.1)
#define     gamma               (0.95)
#define     epsilonm            (1.0)
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

    //float   thpt;
    float n1_s;
    float n2_s;
    float n3_s;
    float self_s;

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

    int                 flag;

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
float           tm, epsilon;
int             debug1, debug2, debug3, debug4, debug5, debug6, sum;
int             n1_pkt, n1_tx, n1_dr, n3_pkt, n2_pkt, n2_tx, n2_dr, n3_tx, n3_dr,n4_pkt, n4_tx, n4_dr, n5_pkt, n5_tx, n5_dr, n6_pkt, n6_tx, n6_dr, n7_pkt, n7_tx, n7_dr, n8_pkt, n8_tx, n8_dr, n9_tx, n9_pkt, n9_dr, n10_pkt, n10_tx, n10_dr, n11_pkt, n11_tx, n11_dr, n12_pkt, n12_tx, n12_dr;
float           th_n1, th_n2, th_n3, rn1, rn2, rn3, rn4, rn5, rn6, rn7, rn8, rn9, rn10, rn11, rn12, act_g, act_g1, act_g2, act_g3, act_g4, act_g5, act_g6, act_g7, act_g8, act_g9, act_g10, act_g11, act_g12;
int             sc_n1, sc_n2, sc_n3, sc_n4, sc_n5, sc_n6, sc_n7, sc_n8, sc_n9, sc_n10, sc_n11, sc_n12;
float           rate1, rate2, rate3, rate4, rate5, rate6, rate7, rate8, rate9, rate10, rate11, rate12, pr1, pr2;
float           rt1, rt2, rt3, rt4, tp1, tp2, c_t, pa1, pa2, pa3, pa4, pa5;
float           true_s1,true_s2,true_s3,true_s4, true_s5,true_s6,true_s7, true_s8, true_s9,true_s10,true_s11, true_s12, true_S, true_st1,true_st2,true_st3,true_st4, true_st5,true_st6,true_st7, true_st8, true_st9,true_st10,true_st11, true_st12;
float           s1,s2,s3,s4, s5, s6, s7, s8, s9, s10, s11, s12, pr_s1, pr_s2, pr_s3, pr_s4, pr_s5, pr_s6, pr_s7, pr_s8, pr_s9, pr_s10, pr_s11, pr_s12, pr_S, s_ratio1, s_ratio2, load_ratio1, load_ratio2;
float           s1_vector[iterations], s2_vector[iterations], s3_vector[iterations],s4_vector[iterations],s5_vector[iterations], s6_vector[iterations],s7_vector[iterations],s8_vector[iterations], s9_vector[iterations],s10_vector[iterations],s11_vector[iterations], s12_vector[iterations], s21_vector[iterations], s23_vector[iterations], s24_vector[iterations], s32_vector[iterations], s34_vector[iterations],  s35_vector[iterations], s42_vector[iterations], s43_vector[iterations], s53_vector[iterations];
int             count_vec[iterations], st1_vector[iterations], st2_vector[iterations], st3_vector[iterations],st4_vector[iterations],st5_vector[iterations];
int             cnt1, cnt2, coll31, coll34, coll42, coll43;
float           fbk1, fbk2, fbk3, pr_fbk1, pr_fbk2, pr_fbk3, d1, d2, d3;
float           prob_in1, prob_in2, prob_in3, prob_in4, prob_in5, prob_in6, prob_in7, prob_in8, prob_in9, prob_in10, prob_in11, prob_in12;
float           est_s12, est_s21, est_s23, est_s24, est_s32, est_s34, est_s35, est_s42, est_s43, est_s53, tpt_31, tpt_34, tpt_42, tpt_43;
float		o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12;
float		pr_o1, pr_o2, pr_o3, pr_o4, pr_o5, pr_o6, pr_o7, pr_o8, pr_o9, pr_o10, pr_o11, pr_o12;
float           s3_from_5, s3_from_4, s2_from_4, s2_from_1, s2_from_3, s4_from_3, s5_from_3, s1_from_2, s3_from_2, s4_from_2, s1_self, s2_self, s3_self, s4_self, s5_self, s6_self, s7_self, s8_self, s9_self, s10_self, s11_self, s12_self;;

/*******************************************************************************/
/* Misc. global variables .......                                              */
/*                                                                             */
/*******************************************************************************/

int             pf_total_gen_packets; /* across the entire network */
int             pf_total_tx_success;
int             new_g;
int             new_g1;
int             new_g2, new_g3, new_g4, new_g5, new_g6, new_g7, new_g8, new_g9, new_g10, new_g11, new_g12;
int             def;
int             coll;
int             def1, coll1, drop1, def3, def4, def5, def6, def7, def8, def9, def10, def11, def12;
int             def2, coll2, drop2, drop3, drop4, drop5, drop6, drop7, drop8, drop9, drop10, drop11, drop12, coll3, coll4, coll5, coll6, coll7, coll8, coll9, coll10, coll11, coll12;
int             collision_count;
int             coll_ct24, coll_ct13, coll_ct34, coll_ct43, coll_ct3, coll_ct4;
int             drop;
int             tx_flag;
int             p_gen1, p_gen2, p_gen3, p_gen4, p_gen5, p_gen6, p_gen7, p_gen8, p_gen9, p_gen10, p_gen11, p_gen12;

int             pkt_ctr1, pkt_ctr2, pkt_ctr3;
int             ctrl, ctrl1;
float           s12, s21, s23, s24, s32, s34, s35, s42, s43, s53;


//float           small_g2;
float           g2[no_of_epochs];
int             window[100];
int             w_id;
float           Q[6][20];
float           Q1[6][20], Q2[6][20], Q3[6][20], Q4[6][20], Q5[6][20], Q6[6][20], Q7[6][20], Q8[6][20], Q9[6][20], Q10[6][20], Q11[6][20];
int             epoch_id;
int             action, action1, action2, action3, action4, action5, action6, action7, action8, action9, action10, action11, action12;
int             s_id;
int             s_id1, s_id2, s_id3, s_id4, s_id5, s_id6, s_id7, s_id8, s_id9, s_id10, s_id11;
int             prev_s_id;
int             prev_s_id1, prev_s_id2, prev_s_id3, prev_s_id4, prev_s_id5, prev_s_id6, prev_s_id7, prev_s_id8, prev_s_id9, prev_s_id10, prev_s_id11;
float           reward;
float           prob_failure;
float           prob_inter;
float           prob_sc;
float           prev_throughput, prev_throughput1;
float           current_throughput;
float           current_throughput1;
float           tmp_var;
int             transmit, transmit1, transmit2, transmit3, transmit4, transmit5, transmit6, transmit7, transmit8, transmit9, transmit10, transmit11;
float           pf1, pf2, pf3, pf4, pf5, pf6, pf7, pf8, pf9, pf10, pf11, pf12, pi1, pi2;
float           delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8, delta9, delta10, delta11, delta12, r1, r2, var1, var2, avg;
float           s_g1, s_g2, s_g3, s_g4, s_g5, s_g6, s_g7, s_g8, s_g9, s_g10, s_g11, s_g12, fct1, fct2,fct3, tn1, tn2, tn3, tn4, tn5, tn;
float           tnn1, tnn42, tnn31, tnn34, tnn43, tnn;
float           fair1, fair2, fair3, fair4, fair5, fair6, fair7, fair8, fair9, fair10, fair11, fair12, pr_fair1, pr_fair2, pr_fair3, pr_fair4, pr_fair5, pr_fair6, pr_fair7, pr_fair8, pr_fair9, pr_fair10, pr_fair11, pr_fair12;
float           fair31, fair32;
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

float           get_exp(float);

int largest(float arr[6][20], int no_of_actions, int state_id);
float maximum(float arr[6][20], int no_of_actions, int state_id);


/*******************************************************************************/
/* function definitions ..                                                     */
/*                                                                             */
/*******************************************************************************/





int largest(float arr[6][20], int no_of_actions, int state_id)
{
    int i;
    int index = 0;
    // Initialize maximum element
    float max = arr[state_id][0];

    // Traverse array elements from second and
    // compare every element with current max
    for (i = 1; i < no_of_actions; i++)
        if (arr[state_id][i] > max){
            max = arr[state_id][i];
            index = i;
        }


    return index;
}


float maximum(float arr[6][20], int no_of_actions, int state_id)
{
    int i;
    int index = 0;
    // Initialize maximum element
    float max = arr[state_id][0];

    // Traverse array elements from second and
    // compare every element with current max
    for (i = 1; i < no_of_actions; i++)
        if (arr[state_id][i] > max){
            max = arr[state_id][i];
            index = i;
        }


    return max;
}






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

        node->n1_s=0.0;
        node->n2_s=0.0;
        node->n3_s=0.0;
        node->self_s=0.0;

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

        event->flag = 0;

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

    /*
    if (time < global_time){
        printf ("%s: .... time < Global time!! Terminating ..", fn);
        Print_Event(event);
        exit (0);
    }
    */

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
    int         k = event->node->Node_ID;
    float      small_g;

    if (k==0){

        n1_pkt++;

        small_g = s_g1;

       /* if (small_g==0){
            return;
        }*/

        event->node->n1_s=s2;

        event->node->self_s=s1_self;

        if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        new_g++;
        new_g1++;
        p_gen1++;

        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }


    float rt1 = rand()/(double)(RAND_MAX);

    if (rt1<=rate1){
        transmit1=1;
    }else{
        transmit1=0;
    }

    if(transmit1==0){
        event->node->dropped_deferred++;
        drop++;
        drop1++;
        n1_dr++;
        new_g=new_g-1;
        new_g1=new_g1-1;
        return;
    }

    /* now deal with the action for this particular event */


    if(transmit1==1){


    n1_tx++;
    /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        //printf("%d\n",event->node->Node_ID);
        def++;
        def1++;
        return;
    }



    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */
    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);
    debug5++;


    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;

        Collision = TRUE;


        int tm1 = Channel->Node_ID;
        if((tm1==1)||(tm1==2)||(tm1==4)||(tm1==6)||(tm1==7)||(tm1==11)){
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

}

    if (k==1){

        n2_pkt++;
        small_g = s_g2;


        event->node->n1_s=s1;
        event->node->n2_s=s3;
        event->node->n3_s=s4;

        event->node->self_s=s2_self;
        //small_g2 = g2[epoch_id];


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen2++;
        new_g++;
        new_g2++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */




    //action= largest(Q,2,s_id);






    float rt2 = rand()/(double)(RAND_MAX);

    if (rt2<=rate2){
        transmit=1;
    }else{
    transmit=0;
    }

    if(transmit==0){
        event->node->dropped_deferred++;
        drop++;
        drop2++;
        debug4++;
        n2_dr++;
        new_g=new_g-1;
        new_g2=new_g2-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit==1){
    debug3++;
    n2_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def2++;
        debug6++;
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);




    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tm2 = Channel->Node_ID;

        if((tm2==0)||(tm2==2)||(tm2==3)||(tm2==4)||(tm2==5)||(tm2==6)||(tm2==7)||(tm2==9)||(tm2==10)){
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
    }





    if (k==2){

        n3_pkt++;
        small_g = small_g3;


        event->node->n1_s=s2;
        event->node->n2_s=s4;
        event->node->n3_s=s5;

        event->node->self_s=s3_self;

        //small_g2 = g2[epoch_id];


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen3++;
        new_g++;
        new_g3++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */




    //action= largest(Q,2,s_id);






    float rt3 = rand()/(double)(RAND_MAX);

    if (rt3<=rate3){
        transmit2=1;
    }else{
    transmit2=0;
    }

    if(transmit2==0){
        event->node->dropped_deferred++;
        drop++;
        drop3++;
        n3_dr++;
        new_g=new_g-1;
        new_g3=new_g3-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit2==1){
    n3_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def3++;
        debug6++;
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);

    /*
    if(ctrl==0){
        pkt_ctr2++;
        if(pkt_ctr2==freq){
            event->flag=1;
            pkt_ctr2=0;
        }
    }
    */


    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        if ((tmk==0)||(tmk==1)||(tmk==3)||(tmk==4)||(tmk==5)||(tmk==7)||(tmk==8)){
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
    }



    if (k==3){

        n4_pkt++;
        small_g = small_g4;


        event->node->n1_s=s2;
        event->node->n2_s=s3;

        event->node->self_s=s4_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen4++;
        new_g++;
        new_g4++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */




    //action= largest(Q,2,s_id);






    float rt4 = rand()/(double)(RAND_MAX);

    if (rt4<=rate4){
        transmit3=1;
    }else{
        transmit3=0;
    }

    if(transmit3==0){
        event->node->dropped_deferred++;
        drop++;
        drop4++;
        n4_dr++;
        new_g=new_g-1;
        new_g4=new_g4-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit3==1){
        n4_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def4++;
        debug6++;
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);

   /*
    if(ctrl1==0){
        pkt_ctr3++;
        if(pkt_ctr3==freq){
            event->flag=1;
            pkt_ctr3=0;
        }
    }
    */


    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

        if ((tmk==1)||(tmk==2)||(tmk==5)||(tmk==7)||(tmk==8)||(tmk==9)||(tmk==11)){
            Collision = TRUE;
            //coll42++;
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
    }




    if (k==4){

        n5_pkt++;
        small_g = small_g5;


        event->node->n1_s=s3;

        event->node->self_s=s5_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen5++;
        new_g++;
        new_g5++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */




    //action= largest(Q,2,s_id);






    float rt5 = rand()/(double)(RAND_MAX);

    if (rt5<=rate5){
        transmit4=1;
    }else{
        transmit4=0;
    }

    if(transmit4==0){
        event->node->dropped_deferred++;
        drop++;
        drop5++;
        n5_dr++;
        new_g=new_g-1;
        new_g5=new_g5-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit4==1){
        n5_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def5++;
        debug6++;
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        if ((tmk==0)||(tmk==1)||(tmk==2)||(tmk==6)||(tmk==7)||(tmk==11)){

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
    }


   if (k==5){

        n6_pkt++;
        small_g = small_g6;


        event->node->n1_s=s3;

        event->node->self_s=s6_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen6++;
        new_g++;
        new_g6++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt6 = rand()/(double)(RAND_MAX);

    if (rt6<=rate6){
        transmit5=1;
    }else{
        transmit5=0;
    }

    if(transmit5==0){
        event->node->dropped_deferred++;
        drop++;
        drop6++;
        n6_dr++;
        new_g=new_g-1;
        new_g6=new_g6-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit5==1){
        n6_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def6++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

        if ((tmk==1)||(tmk==2)||(tmk==3)||(tmk==6)||(tmk==7)||(tmk==8)||(tmk==9)||(tmk==10)){

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
    }

   if (k==6){

        n7_pkt++;
        small_g = small_g7;


        event->node->n1_s=s3;

        event->node->self_s=s7_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen7++;
        new_g++;
        new_g7++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt7 = rand()/(double)(RAND_MAX);

    if (rt7<=rate7){
        transmit6=1;
    }else{
        transmit6=0;
    }

    if(transmit6==0){
        event->node->dropped_deferred++;
        drop++;
        drop7++;
        n7_dr++;
        new_g=new_g-1;
        new_g7=new_g7-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit6==1){
        n7_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def7++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

        if ((tmk==0)||(tmk==1)||(tmk==4)||(tmk==5)||(tmk==7)||(tmk==8)||(tmk==9)||(tmk==11)||(tmk==10)){
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
    }

   if (k==7){

        n8_pkt++;
        small_g = small_g8;


        event->node->n1_s=s3;

        event->node->self_s=s8_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen8++;
        new_g++;
        new_g8++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt8 = rand()/(double)(RAND_MAX);

    if (rt8<=rate8){
        transmit7=1;
    }else{
        transmit7=0;
    }

    if(transmit7==0){
        event->node->dropped_deferred++;
        drop++;
        drop8++;
        n8_dr++;
        new_g=new_g-1;
        new_g8=new_g8-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit7==1){
        n8_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def8++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        Collision = TRUE;


	/*
        if ((tmk==0)||(tmk==1)||(tmk==2)||(tmk==3)||(tmk==4)||(tmk==5)||(tmk==6)||(tmk==8)||(tmk==9)||(tmk==10)){
            Collision = TRUE;
        }
	*/


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
    }

    if (k==8){

        n9_pkt++;
        small_g = small_g9;


        event->node->n1_s=s3;

        event->node->self_s=s9_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen9++;
        new_g++;
        new_g9++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt9 = rand()/(double)(RAND_MAX);

    if (rt9<=rate9){
        transmit8=1;
    }else{
        transmit8=0;
    }

    if(transmit8==0){
        event->node->dropped_deferred++;
        drop++;
        drop9++;
        n9_dr++;
        new_g=new_g-1;
        new_g9=new_g9-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit8==1){
        n9_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def9++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

        if ((tmk==2)||(tmk==3)||(tmk==5)||(tmk==6)||(tmk==7)||(tmk==9)||(tmk==10)||(tmk==11)){
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
    }


if (k==9){

        n10_pkt++;
        small_g = small_g10;


        event->node->n1_s=s3;

        event->node->self_s=s10_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen10++;
        new_g++;
        new_g10++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt10 = rand()/(double)(RAND_MAX);

    if (rt10<=rate10){
        transmit9=1;
    }else{
        transmit9=0;
    }

    if(transmit9==0){
        event->node->dropped_deferred++;
        drop++;
        drop10++;
        n10_dr++;
        new_g=new_g-1;
        new_g10=new_g10-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit9==1){
        n10_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def10++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

         if ((tmk==1)||(tmk==3)||(tmk==5)||(tmk==6)||(tmk==7)||(tmk==8)||(tmk==10)||(tmk==11)){
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
    }


if (k==10){

        n11_pkt++;
        small_g = small_g11;


        event->node->n1_s=s3;

        event->node->self_s=s11_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen11++;
        new_g++;
        new_g11++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt11 = rand()/(double)(RAND_MAX);

    if (rt11<=rate11){
        transmit10=1;
    }else{
        transmit10=0;
    }

    if(transmit10==0){
        event->node->dropped_deferred++;
        drop++;
        drop11++;
        n11_dr++;
        new_g=new_g-1;
        new_g11=new_g11-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit10==1){
        n11_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def11++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

         if ((tmk==1)||(tmk==5)||(tmk==4)||(tmk==0)||(tmk==7)||(tmk==9)||(tmk==6)||(tmk==8)||(tmk==11)){
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
    }


if (k==11){

        n12_pkt++;
        small_g = small_g12;


        event->node->n1_s=s3;

        event->node->self_s=s12_self;


    if (event->Event_type == Ev_P_GEN){

        /* this is a generation event; not deferred transmission */
        /* among other things, schedule another generation       */

        event->node->generated++;

        pf_total_gen_packets++;
        p_gen12++;
        new_g++;
        new_g12++;





        /* schedule the next generation event for this node */
        Insert_Event(Ev_P_GEN, event->Event_time + get_exp(1.0/((float)small_g/(double)packet_duration)), event->node);
    }

    /* now deal with the action for this particular event */







    float rt12 = rand()/(double)(RAND_MAX);

    if (rt12<=rate12){
        transmit11=1;
    }else{
        transmit11=0;
    }

    if(transmit11==0){
        event->node->dropped_deferred++;
        drop++;
        drop12++;
        n12_dr++;
        new_g=new_g-1;
        new_g12=new_g12-1;
        return;
    }





    /* node is not in Tx; ready to make a Tx for the current packet based on P_WS etc. */

    /* schedule an end-of-Tx event */




    if(transmit11==1){
        n12_tx++;


        /* if the node is currently within a Tx, drop/defer the packet */
    if (event->node->Tx_Status == TRUE){

        event->node->dropped_deferred++;
        def++;
        def12++;
        debug6++; //....
        return;
    }

    Insert_Event(Ev_End_Tx, event->Event_time + packet_duration, event->node);



    /* insert this node in the beginning of the channel list, which contains a list of nodes on the channel  */
    if (!Channel){  // channel is free right now
        event->node->Next = (NodeInfo *)NULL;
        event->node->Prev = (NodeInfo *)NULL;
    }
    else{ // channel has one or more nodes already ..
        event->node->Next = Channel;
        Channel->Prev = event->node;
        int tmk = Channel->Node_ID;
        //Collision = TRUE;

         if ((tmk==0)||(tmk==3)||(tmk==4)||(tmk==6)||(tmk==8)||(tmk==7)||(tmk==9)||(tmk==10)){
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
        }
    }










}

/*  handler/action routine for end of packet transmission */
void    End_Transmission_Action(Event *event)
{
    if (!event->node->Coll_Status){
        event->node->success++;
        pf_total_tx_success++;


        if(event->node->Node_ID==0){
            sc_n1++;
            s2_from_1= event->node->n1_s;
            //s12 = event->node->self_s;
        }
        if(event->node->Node_ID==1){
            sc_n2++;
            s1_from_2= event->node->n1_s;
            s3_from_2= event->node->n2_s;
            s4_from_2= event->node->n3_s;
            s21 = event->node->self_s;
            s23 = event->node->self_s;
            s24 = event->node->self_s;
        }
        if(event->node->Node_ID==2){
            sc_n3++;
            s2_from_3= event->node->n1_s;
            s4_from_3= event->node->n2_s;
            s5_from_3= event->node->n3_s;
            s32 = event->node->self_s;
            s34 = event->node->self_s;
            s35 = event->node->self_s;
        }

        if(event->node->Node_ID==3){
            sc_n4++;
            s2_from_4= event->node->n1_s;
            s3_from_4= event->node->n2_s;
            s42 = event->node->self_s;
            s43 = event->node->self_s;
        }
        if(event->node->Node_ID==4){
            sc_n5++;
            s3_from_5= event->node->n1_s;
            s53 = event->node->self_s;
        }
        if(event->node->Node_ID==5){
            sc_n6++;
        }
        if(event->node->Node_ID==6){
            sc_n7++;
        }
        if(event->node->Node_ID==7){
            sc_n8++;
        }
        if(event->node->Node_ID==8){
            sc_n9++;
        }
        if(event->node->Node_ID==9){
            sc_n10++;
        }
        if(event->node->Node_ID==10){
            sc_n11++;
        }
        if(event->node->Node_ID==11){
            sc_n12++;
        }


    }else{
        event->node->collided++;
        coll++;
        if(event->node->Node_ID==0){
            coll1++;
        }
        if(event->node->Node_ID==1){
            coll2++;
        }
        if(event->node->Node_ID==2){
            coll3++;
        }
        if(event->node->Node_ID==3){
            coll4++;
        }
        if(event->node->Node_ID==4){
            coll5++;
        }
        if(event->node->Node_ID==5){
            coll6++;
        }
        if(event->node->Node_ID==6){
            coll7++;
        }
        if(event->node->Node_ID==7){
            coll8++;
        }
        if(event->node->Node_ID==8){
            coll9++;
        }
        if(event->node->Node_ID==9){
            coll10++;
        }
        if(event->node->Node_ID==10){
            coll11++;
        }
        if(event->node->Node_ID==11){
            coll12++;
        }
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
    NodeInfo        *node1;
    int             l;
    //printf("Epoch ID: %d\n", epoch_id);








    while (!Sim_Stop){

        event = Remove_Front_Event();
        global_time = event->Event_time;

        //printf("Epoch ID: %d\n", epoch_id);

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

        if ((global_time >= Simulation_Duration)||(pf_total_gen_packets>packets_per_epoch))
            Sim_Stop = TRUE;


//            printf("Global Time : %d\n",global_time);



    }



    prob_failure = (double)(def)/(double)(pf_total_gen_packets-drop);
    pf1 = (double)(def1)/(double)(p_gen1-drop1+0.000001);
    pf2 = (double)(def2)/(double)(p_gen2-drop2+0.000001);
    pf3 = (double)(def3)/(double)(p_gen3-drop3+0.000001);
    pf4 = (double)(def4)/(double)(p_gen4-drop4+0.000001);
    pf5 = (double)(def5)/(double)(p_gen5-drop5+0.000001);
    pf6 = (double)(def6)/(double)(p_gen6-drop6+0.000001);
    pf7 = (double)(def7)/(double)(p_gen7-drop7+0.000001);
    pf8 = (double)(def8)/(double)(p_gen8-drop8+0.000001);
    pf9 = (double)(def9)/(double)(p_gen9-drop9+0.000001);
    pf10 = (double)(def10)/(double)(p_gen10-drop10+0.000001);
    pf11 = (double)(def11)/(double)(p_gen11-drop11+0.000001);
    pf12 = (double)(def12)/(double)(p_gen12-drop12+0.000001);
    //prob_failure = (100-(float)sum)/100;

    prob_inter = (double)(coll)/(double)(pf_total_gen_packets-drop+0.000001);
    prob_in1 = (double)(coll1)/(double)(p_gen1-drop1+0.000001);
    prob_in2 = (double)(coll2)/(double)(p_gen2-drop2+0.000001);
    prob_in3 = (double)(coll3)/(double)(p_gen3-drop3+0.000001);
    prob_in4 = (double)(coll4)/(double)(p_gen4-drop4+0.000001);
    prob_in5 = (double)(coll5)/(double)(p_gen5-drop5+0.000001);
    prob_in6 = (double)(coll6)/(double)(p_gen6-drop6+0.000001);
    prob_in7 = (double)(coll7)/(double)(p_gen7-drop7+0.000001);
    prob_in8 = (double)(coll8)/(double)(p_gen8-drop8+0.000001);
    prob_in9 = (double)(coll9)/(double)(p_gen9-drop9+0.000001);
    prob_in10 = (double)(coll10)/(double)(p_gen10-drop10+0.000001);
    prob_in11 = (double)(coll11)/(double)(p_gen11-drop11+0.000001);
    prob_in12 = (double)(coll12)/(double)(p_gen12-drop12+0.000001);
    //pi1 = (double)(def1)/(double)(p_gen1-drop1);

    prev_s_id = s_id;
    prev_s_id1 = s_id1;
    prev_s_id2 = s_id2;
    prev_s_id3 = s_id3;




    if ((new_g==0)){
        current_throughput = ((1.00-pkt_er_rate)*(float)pf_total_tx_success/(float)(pf_total_gen_packets));
    }if ((new_g!=0)){
        current_throughput = ((1.00-pkt_er_rate)*(float)pf_total_tx_success/(float)(new_g));
    }


    current_throughput1 = ((1.00-pkt_er_rate)*(float)pf_total_tx_success/(float)(pf_total_gen_packets));



    th_n1 = (float)sc_n1/(float)pf_total_gen_packets;
    th_n2 = (float)sc_n2/(float)pf_total_gen_packets;
    th_n3 = (float)sc_n3/(float)pf_total_gen_packets;


    pa1 = act_g1;
    pa2 = act_g2;
    pa3 = act_g3;

    act_g = ((double)(new_g)*(double)packet_duration)/(float)global_time;
    act_g1 = ((double)(new_g1)*(double)packet_duration)/(float)global_time;
    act_g2 = ((double)(new_g2)*(double)packet_duration)/(float)global_time;
    act_g3 = ((double)(new_g3)*(double)packet_duration)/(float)global_time;
    act_g4 = ((double)(new_g4)*(double)packet_duration)/(float)global_time;
    act_g5 = ((double)(new_g5)*(double)packet_duration)/(float)global_time;
    act_g6 = ((double)(new_g6)*(double)packet_duration)/(float)global_time;
    act_g7 = ((double)(new_g7)*(double)packet_duration)/(float)global_time;
    act_g8 = ((double)(new_g8)*(double)packet_duration)/(float)global_time;
    act_g9 = ((double)(new_g9)*(double)packet_duration)/(float)global_time;
    act_g10 = ((double)(new_g10)*(double)packet_duration)/(float)global_time;
    act_g11 = ((double)(new_g11)*(double)packet_duration)/(float)global_time;
    act_g12 = ((double)(new_g12)*(double)packet_duration)/(float)global_time;




    tn = (act_g)*(1.00-pkt_er_rate)*(double)pf_total_tx_success/(double)(new_g+0.000001);
    s1= (act_g1)*(1.00-pkt_er_rate)*(double)sc_n1/(double)(new_g1+0.000001);
    s2= (act_g2)*(1.00-pkt_er_rate)*(double)sc_n2/(double)(new_g2+0.000001);
    s3= (act_g3)*(1.00-pkt_er_rate)*(double)sc_n3/(double)(new_g3+0.000001);
    s4= (act_g4)*(1.00-pkt_er_rate)*(double)sc_n4/(double)(new_g4+0.000001);
    s5= (act_g5)*(1.00-pkt_er_rate)*(double)sc_n5/(double)(new_g5+0.000001);
    s6= (act_g6)*(1.00-pkt_er_rate)*(double)sc_n6/(double)(new_g6+0.000001);
    s7= (act_g7)*(1.00-pkt_er_rate)*(double)sc_n7/(double)(new_g7+0.000001);
    s8= (act_g8)*(1.00-pkt_er_rate)*(double)sc_n8/(double)(new_g8+0.000001);
    s9= (act_g9)*(1.00-pkt_er_rate)*(double)sc_n9/(double)(new_g9+0.000001);
    s10= (act_g10)*(1.00-pkt_er_rate)*(double)sc_n10/(double)(new_g10+0.000001);
    s11= (act_g11)*(1.00-pkt_er_rate)*(double)sc_n11/(double)(new_g11+0.000001);
    s12= (act_g12)*(1.00-pkt_er_rate)*(double)sc_n12/(double)(new_g12+0.000001);


    tpt_31=(act_g3)*((1.00-pkt_er_rate)*(double)sc_n3+(double)coll34)/(double)(new_g3+0.000001);
    tpt_34=(act_g3)*((1.00-pkt_er_rate)*(double)sc_n3+(double)coll31)/(double)(new_g3+0.000001);
    tpt_42=(act_g4)*((1.00-pkt_er_rate)*(double)sc_n4+(double)coll43)/(double)(new_g4+0.000001);
    tpt_43=(act_g4)*((1.00-pkt_er_rate)*(double)sc_n4+(double)coll42)/(double)(new_g4+0.000001);




    s_ratio1=s1/s2;
    s_ratio2=s3/s4;

    //tnn = (act_g)*(double)pf_total_tx_success/(double)(new_g+0.000001);
    //tnn1= (act_g1)*(double)(sc_n1+coll_ct3)/(double)(new_g1+0.000001);
    tnn31= (act_g3)*(1.00-pkt_er_rate)*(double)(sc_n3+coll_ct34)/(double)(new_g3+0.000001);
    tnn42= (act_g4)*(1.00-pkt_er_rate)*(double)(sc_n4+coll_ct34)/(double)(new_g4+0.000001);
    tnn43= (act_g4)*(1.00-pkt_er_rate)*(double)(sc_n4+coll_ct24)/(double)(new_g4+0.000001);
    tnn34= (act_g3)*(1.00-pkt_er_rate)*(double)(sc_n3+coll_ct13)/(double)(new_g3+0.000001);



        prev_throughput1=current_throughput1;
        prev_throughput = current_throughput;




    if (prob_in1<0.165){
            s_id = 0;
        }
        if ((prob_in1>=0.165)&&(prob_in1<0.33)){
            s_id = 1;
        }
        if ((prob_in1>=0.33)&&(prob_in1<0.495)){
            s_id = 2;
        }
        if ((prob_in1>=0.495)&&(prob_in1<0.66)){
            s_id = 3;
        }
        if ((prob_in1>=0.66)&&(prob_in1<0.825)){
            s_id = 4;
        }
        if (prob_in1>=0.825){
            s_id = 5;
        }

        if (prob_in2<0.165){
            s_id1 = 0;
        }
        if ((prob_in2>=0.165)&&(prob_in2<0.33)){
            s_id1 = 1;
        }
        if ((prob_in2>=0.33)&&(prob_in2<0.495)){
            s_id1 = 2;
        }
        if ((prob_in2>=0.495)&&(prob_in2<0.66)){
            s_id1 = 3;
        }
        if ((prob_in2>=0.66)&&(prob_in2<0.825)){
            s_id1 = 4;
        }
        if (prob_in2>=0.825){
            s_id1 = 5;
        }


        if (prob_in3<0.165){
            s_id2 = 0;
        }
        if ((prob_in3>=0.165)&&(prob_in3<0.33)){
            s_id2 = 1;
        }
        if ((prob_in3>=0.33)&&(prob_in3<0.495)){
            s_id2 = 2;
        }
        if ((prob_in3>=0.495)&&(prob_in3<0.66)){
            s_id2 = 3;
        }
        if ((prob_in3>=0.66)&&(prob_in3<0.825)){
            s_id2 = 4;
        }
        if (prob_in3>=0.825){
            s_id2 = 5;
        }

        if (prob_in4<0.165){
            s_id3 = 0;
        }
        if ((prob_in4>=0.165)&&(prob_in4<0.33)){
            s_id3 = 1;
        }
        if ((prob_in4>=0.33)&&(prob_in4<0.495)){
            s_id3 = 2;
        }
        if ((prob_in4>=0.495)&&(prob_in4<0.66)){
            s_id3 = 3;
        }
        if ((prob_in4>=0.66)&&(prob_in4<0.825)){
            s_id3 = 4;
        }
        if (prob_in4>=0.825){
            s_id3 = 5;
        }

        if (prob_in5<0.165){
            s_id4 = 0;
        }
        if ((prob_in5>=0.165)&&(prob_in5<0.33)){
            s_id4 = 1;
        }
        if ((prob_in5>=0.33)&&(prob_in5<0.495)){
            s_id4 = 2;
        }
        if ((prob_in5>=0.495)&&(prob_in5<0.66)){
            s_id4 = 3;
        }
        if ((prob_in5>=0.66)&&(prob_in5<0.825)){
            s_id4 = 4;
        }
        if (prob_in5>=0.825){
            s_id4 = 5;
        }

        if (prob_in6<0.165){
            s_id5 = 0;
        }
        if ((prob_in6>=0.165)&&(prob_in6<0.33)){
            s_id5 = 1;
        }
        if ((prob_in6>=0.33)&&(prob_in6<0.495)){
            s_id5 = 2;
        }
        if ((prob_in6>=0.495)&&(prob_in6<0.66)){
            s_id5 = 3;
        }
        if ((prob_in6>=0.66)&&(prob_in6<0.825)){
            s_id5 = 4;
        }
        if (prob_in6>=0.825){
            s_id5 = 5;
        }

        if (prob_in7<0.165){
            s_id6 = 0;
        }
        if ((prob_in7>=0.165)&&(prob_in7<0.33)){
            s_id6 = 1;
        }
        if ((prob_in7>=0.33)&&(prob_in7<0.495)){
            s_id6 = 2;
        }
        if ((prob_in7>=0.495)&&(prob_in7<0.66)){
            s_id6 = 3;
        }
        if ((prob_in7>=0.66)&&(prob_in7<0.825)){
            s_id6 = 4;
        }
        if (prob_in7>=0.825){
            s_id6 = 5;
        }

        if (prob_in8<0.165){
            s_id7 = 0;
        }
        if ((prob_in8>=0.165)&&(prob_in8<0.33)){
            s_id7 = 1;
        }
        if ((prob_in8>=0.33)&&(prob_in8<0.495)){
            s_id7 = 2;
        }
        if ((prob_in8>=0.495)&&(prob_in8<0.66)){
            s_id7 = 3;
        }
        if ((prob_in8>=0.66)&&(prob_in8<0.825)){
            s_id7 = 4;
        }
        if (prob_in8>=0.825){
            s_id7 = 5;
        }

        if (prob_in9<0.165){
            s_id8 = 0;
        }
        if ((prob_in9>=0.165)&&(prob_in9<0.33)){
            s_id8 = 1;
        }
        if ((prob_in9>=0.33)&&(prob_in9<0.495)){
            s_id8 = 2;
        }
        if ((prob_in9>=0.495)&&(prob_in9<0.66)){
            s_id8 = 3;
        }
        if ((prob_in9>=0.66)&&(prob_in9<0.825)){
            s_id8 = 4;
        }
        if (prob_in9>=0.825){
            s_id8 = 5;
        }


        if (prob_in10<0.165){
            s_id9 = 0;
        }
        if ((prob_in10>=0.165)&&(prob_in10<0.33)){
            s_id9 = 1;
        }
        if ((prob_in10>=0.33)&&(prob_in10<0.495)){
            s_id9 = 2;
        }
        if ((prob_in10>=0.495)&&(prob_in10<0.66)){
            s_id9 = 3;
        }
        if ((prob_in10>=0.66)&&(prob_in10<0.825)){
            s_id9 = 4;
        }
        if (prob_in10>=0.825){
            s_id9 = 5;
        }

        if (prob_in11<0.165){
            s_id10 = 0;
        }
        if ((prob_in11>=0.165)&&(prob_in11<0.33)){
            s_id10 = 1;
        }
        if ((prob_in11>=0.33)&&(prob_in11<0.495)){
            s_id10 = 2;
        }
        if ((prob_in11>=0.495)&&(prob_in11<0.66)){
            s_id10 = 3;
        }
        if ((prob_in11>=0.66)&&(prob_in11<0.825)){
            s_id10 = 4;
        }
        if (prob_in11>=0.825){
            s_id10 = 5;
        }

        if (prob_in12<0.165){
            s_id11 = 0;
        }
        if ((prob_in12>=0.165)&&(prob_in12<0.33)){
            s_id11 = 1;
        }
        if ((prob_in12>=0.33)&&(prob_in12<0.495)){
            s_id11 = 2;
        }
        if ((prob_in12>=0.495)&&(prob_in12<0.66)){
            s_id11 = 3;
        }
        if ((prob_in12>=0.66)&&(prob_in12<0.825)){
            s_id11 = 4;
        }
        if (prob_in12>=0.825){
            s_id11 = 5;
        }












    if (act_g1>s_g1){
        act_g1=s_g1;
    }

    if (act_g2>s_g2){
        act_g2=s_g2;
    }
    if (act_g3>s_g3){
        act_g3=s_g3;
    }
    if (act_g4>s_g4){
        act_g4=s_g4;
    }
    if (act_g5>s_g5){
        act_g5=s_g5;
    }
    if (act_g6>s_g6){
        act_g6=s_g6;
    }
    if (act_g7>s_g7){
        act_g7=s_g7;
    }
    if (act_g8>s_g8){
        act_g8=s_g8;
    }
    if (act_g9>s_g9){
        act_g9=s_g9;
    }
    if (act_g10>s_g10){
        act_g10=s_g10;
    }
    if (act_g11>s_g11){
        act_g11=s_g11;
    }
    if (act_g12>s_g12){
        act_g12=s_g12;
    }


    if (act_g>s_g2+s_g1+s_g3+s_g4+s_g5+s_g6+s_g7+s_g8+s_g9+s_g10+s_g11+s_g12){
        act_g=s_g2+s_g1+s_g3+s_g4+s_g5+s_g6+s_g7+s_g8+s_g9+s_g10+s_g11+s_g12;
    }

    tp1=((1.00-pkt_er_rate)*(double)pf_total_tx_success/(double)global_time)*(double)packet_duration;
    tp2 = (act_g)*(1.00-pkt_er_rate)*(double)pf_total_tx_success/(double)(new_g+0.000001);



}

/* Main function that initializes things and starts simulation and then prints stats .. */

double randZeroToOne()
{
    return (1)*rand()/(double)(RAND_MAX);
}


int main()
{
    char    *fn="main";
    int      iter, iter1, iter2;
    int      i;
    FILE    *fpt, *fpt1;

    fpt = fopen("rev1_nfc_12_n_1_test_0_3_0_1_c.txt","w");
    fpt1 = fopen("rev1_intm_12_n_14_testh_0_1.txt","w");

    srand(time(0));
#ifdef S_DEBUG
    printf ("I am here ... \n");
#endif


    for(i=0;i<no_of_epochs;i++){
        g2[i]=10*(float)i/(float)(no_of_epochs);
        //g2[1]=0;
    }



        for(iter=0;iter<6;iter++){
        for(iter2=0;iter2<20;iter2++){
            Q[iter][iter2] = 1.0;
            Q1[iter][iter2] = 1.0;
            Q2[iter][iter2] = 1.0;
            Q3[iter][iter2] = 1.0;
            Q4[iter][iter2] = 1.0;
            Q5[iter][iter2] = 1.0;
            Q6[iter][iter2] = 1.0;
            Q7[iter][iter2] = 1.0;
            Q8[iter][iter2] = 1.0;
            Q9[iter][iter2] = 1.0;
            Q10[iter][iter2] = 1.0;
            Q11[iter][iter2] = 1.0;
        }
    }




    rate1=1.0;
    rate2=1.0;
    rate3=1.0;
    rate4=1.0;
    rate5=1.0;
    rate6=1.0;
    rate7=1.0;
    rate8=1.0;
    rate9=1.0;
    rate10=1.0;
    rate11=1.0;
    rate12=1.0;
    //rate4=1.0;

    pr1=0.0;
    pr2=0.0;

    s_g1=small_g1;
    s_g2=small_g2;
    s_g3=small_g3;
    s_g4=small_g4;
    s_g5=small_g5;
    s_g6=small_g6;
    s_g7=small_g7;
    s_g8=small_g8;
    s_g9=small_g9;
    s_g10=small_g10;
    s_g11=small_g11;
    s_g12=small_g12;
    //s_g4=small_g4;




    Create_Nodes();
    Create_Ev_Pool();
    w_id = 0;
    prob_failure = 0.0;
    prob_inter = 0.0;
    prob_in1 = 0.0;
    prob_in2 = 0.0;
    prob_in3 = 0.0;
    prob_in4 = 0.0;
    prob_in5 = 0.0;
    prob_in6 = 0.0;
    prob_in7 = 0.0;
    prob_in8 = 0.0;
    prob_in9 = 0.0;
    prob_in10 = 0.0;
    prob_in11 = 0.0;
    prob_in12 = 0.0;
    s_id=0;
    s_id1=0;
    s_id2=0;
    s_id3=0;
    s_id4=0;
    s_id5=0;
    s_id6=0;
    s_id7=0;
    s_id8=0;
    s_id9=0;
    s_id10=0;
    s_id11=0;
    s1=0.0;
    s2=0.0;
    s3=0.0;
    s4=0.0;
    s5=0.0;
    s6=0.0;
    s7=0.0;
    s8=0.0;
    s9=0.0;
    s10=0.0;
    s11=0.0;
    s12=0.0;
    prev_throughput=0;
    prev_throughput1=0;
    current_throughput=0;
    current_throughput1=0;
    th_n1= 0.0;
    th_n2 = 0.0;
    th_n3 = 0.0;
    tmp_var=0;
    cnt1=0;
    cnt2=0;
    pkt_ctr2=0;
    //s31=0.0;
    //s34=0.0;
    //s42=0.0;
    //s43=0.0;
    //srand(time(NULL));
    for(epoch_id=0;epoch_id<no_of_epochs;epoch_id++)
    {



        //ctrl=0;
        //ctrl1=0;
        //pkt_ctr2=0;
        //pkt_ctr3=0;
        //printf("%d",ctrl);


        prev_s_id=true_st1;
        prev_s_id1=true_st2;
        prev_s_id2=true_st3;
        prev_s_id3=true_st4;
        prev_s_id4=true_st5;
        prev_s_id5=true_st6;
        prev_s_id6=true_st7;
        prev_s_id7=true_st8;
        prev_s_id8=true_st9;
        prev_s_id9=true_st10;
        prev_s_id10=true_st11;
        prev_s_id11=true_st12;


    epsilon = (float)epsilonm*(float)exp(-((float)epoch_id/2200.0));

    float random = rand()/(double)(RAND_MAX);
    float random1 = rand()/(double)(RAND_MAX);
    float random2 = rand()/(double)(RAND_MAX);


    if (random<epsilon){
        if (random1<0.05){
            action1 = 0;
        }
        if ((random1>=0.05)&&(random1<0.1)){
            action1 = 1;
        }
        if ((random1>=0.10)&&(random1<0.15)){
            action1 = 2;
        }
        if ((random1>=0.15)&&(random1<0.2)){
            action1 = 3;
        }
        if ((random1>=0.2)&&(random1<0.25)){
            action1 = 4;
        }
        if ((random1>=0.25)&&(random1<0.3)){
            action1 = 5;
        }
        if ((random1>=0.3)&&(random1<0.35)){
            action1 = 6;
        }
        if ((random1>=0.35)&&(random1<0.4)){
            action1 = 7;
        }
        if ((random1>=0.4)&&(random1<0.45)){
            action1 = 8;
        }
        if ((random1>=0.45)&&(random1<0.5)){
            action1 = 9;
        }
        if ((random1>=0.5)&&(random1<0.55)){
            action1 = 10;
        }
        if ((random1>=0.55)&&(random1<0.6)){
            action1 = 11;
        }
        if ((random1>=0.6)&&(random1<0.65)){
            action1 = 12;
        }
        if ((random1>=0.65)&&(random1<0.7)){
            action1 = 13;
        }
        if ((random1>=0.7)&&(random1<0.75)){
            action1 = 14;
        }
        if ((random1>=0.75)&&(random1<0.8)){
            action1 = 15;
        }
        if ((random1>=0.8)&&(random1<0.85)){
            action1 = 16;
        }
        if ((random1>=0.85)&&(random1<0.9)){
            action1 = 17;
        }
        if ((random1>=0.9)&&(random1<0.95)){
            action1 = 18;
        }
        if ((random1>=0.95)){
            action1 = 18;
        }
    }
    else{
        action1 = largest(Q,20,true_st1);
    }


    if(action1==0){
       rate1=0.05;
    }
    if(action1==1){
        rate1=0.1;
    }
    if(action1==2){
        rate1=0.15;
    }
    if(action1==3){
       rate1=0.2;
    }
    if(action1==4){
        rate1=0.25;
    }
    if(action1==5){
        rate1=0.3;
    }
    if(action1==6){
       rate1=0.35;
    }
    if(action1==7){
        rate1=0.4;
    }
    if(action1==8){
        rate1=0.45;
    }
    if(action1==9){
        rate1=0.5;
    }
    if(action1==10){
       rate1=0.55;
    }
    if(action1==11){
        rate1=0.6;
    }
    if(action1==12){
        rate1=0.65;
    }
    if(action1==13){
       rate1=0.7;
    }
    if(action1==14){
        rate1=0.75;
    }
    if(action1==15){
        rate1=0.8;
    }
    if(action1==16){
       rate1=0.85;
    }
    if(action1==17){
        rate1=0.9;
    }
    if(action1==18){
        rate1=0.95;
    }
    if(action1==19){
        rate1=1.00;
    }



    pr2=rate2;




    float randoMA = rand()/(double)(RAND_MAX);
    float randoMA1 = rand()/(double)(RAND_MAX);
    float randoMA2 = rand()/(double)(RAND_MAX);


    if (randoMA<epsilon){
        if (randoMA1<0.05){
            action2 = 0;
        }
        if ((randoMA1>=0.05)&&(randoMA1<0.1)){
            action2 = 1;
        }
        if ((randoMA1>=0.10)&&(randoMA1<0.15)){
            action2 = 2;
        }
        if ((randoMA1>=0.15)&&(randoMA1<0.2)){
            action2 = 3;
        }
        if ((randoMA1>=0.2)&&(randoMA1<0.25)){
            action2 = 4;
        }
        if ((randoMA1>=0.25)&&(randoMA1<0.3)){
            action2 = 5;
        }
        if ((randoMA1>=0.3)&&(randoMA1<0.35)){
            action2 = 6;
        }
        if ((randoMA1>=0.35)&&(randoMA1<0.4)){
            action2 = 7;
        }
        if ((randoMA1>=0.4)&&(randoMA1<0.45)){
            action2 = 8;
        }
        if ((randoMA1>=0.45)&&(randoMA1<0.5)){
            action2 = 9;
        }
        if ((randoMA1>=0.5)&&(randoMA1<0.55)){
            action2 = 10;
        }
        if ((randoMA1>=0.55)&&(randoMA1<0.6)){
            action2 = 11;
        }
        if ((randoMA1>=0.6)&&(randoMA1<0.65)){
            action2 = 12;
        }
        if ((randoMA1>=0.65)&&(randoMA1<0.7)){
            action2 = 13;
        }
        if ((randoMA1>=0.7)&&(randoMA1<0.75)){
            action2 = 14;
        }
        if ((randoMA1>=0.75)&&(randoMA1<0.8)){
            action2 = 15;
        }
        if ((randoMA1>=0.8)&&(randoMA1<0.85)){
            action2 = 16;
        }
        if ((randoMA1>=0.85)&&(randoMA1<0.9)){
            action2 = 17;
        }
        if ((randoMA1>=0.9)&&(randoMA1<0.95)){
            action2 = 18;
        }
        if ((randoMA1>=0.95)){
            action2 = 19;
        }
    }
    else{
        action2 = largest(Q1,20,true_st2);
    }


    if(action2==0){
       rate2=0.05;
    }
    if(action2==1){
        rate2=0.1;
    }
    if(action2==2){
        rate2=0.15;
    }
    if(action2==3){
       rate2=0.2;
    }
    if(action2==4){
        rate2=0.25;
    }
    if(action2==5){
        rate2=0.3;
    }
    if(action2==6){
       rate2=0.35;
    }
    if(action2==7){
        rate2=0.4;
    }
    if(action2==8){
        rate2=0.45;
    }
    if(action2==9){
        rate2=0.5;
    }
    if(action2==10){
       rate2=0.55;
    }
    if(action2==11){
        rate2=0.6;
    }
    if(action2==12){
        rate2=0.65;
    }
    if(action2==13){
       rate2=0.7;
    }
    if(action2==14){
        rate2=0.75;
    }
    if(action2==15){
        rate2=0.8;
    }
    if(action2==16){
       rate2=0.85;
    }
    if(action2==17){
        rate2=0.9;
    }
    if(action2==18){
        rate2=0.95;
    }
    if(action2==19){
        rate2=1.00;
    }





    float randomb = rand()/(double)(RAND_MAX);
    float randomb1 = rand()/(double)(RAND_MAX);
    float randomb2 = rand()/(double)(RAND_MAX);


    if (randomb<epsilon){
        if (randomb1<0.05){
            action3 = 0;
        }
        if ((randomb1>=0.05)&&(randomb1<0.1)){
            action3 = 1;
        }
        if ((randomb1>=0.10)&&(randomb1<0.15)){
            action3 = 2;
        }
        if ((randomb1>=0.15)&&(randomb1<0.2)){
            action3 = 3;
        }
        if ((randomb1>=0.2)&&(randomb1<0.25)){
            action3 = 4;
        }
        if ((randomb1>=0.25)&&(randomb1<0.3)){
            action3 = 5;
        }
        if ((randomb1>=0.3)&&(randomb1<0.35)){
            action3 = 6;
        }
        if ((randomb1>=0.35)&&(randomb1<0.4)){
            action3 = 7;
        }
        if ((randomb1>=0.4)&&(randomb1<0.45)){
            action3 = 8;
        }
        if ((randomb1>=0.45)&&(randomb1<0.5)){
            action3 = 9;
        }
        if ((randomb1>=0.5)&&(randomb1<0.55)){
            action3 = 10;
        }
        if ((randomb1>=0.55)&&(randomb1<0.6)){
            action3 = 11;
        }
        if ((randomb1>=0.6)&&(randomb1<0.65)){
            action3 = 12;
        }
        if ((randomb1>=0.65)&&(randomb1<0.7)){
            action3 = 13;
        }
        if ((randomb1>=0.7)&&(randomb1<0.75)){
            action3 = 14;
        }
        if ((randomb1>=0.75)&&(randomb1<0.8)){
            action3 = 15;
        }
        if ((randomb1>=0.8)&&(randomb1<0.85)){
            action3 = 16;
        }
        if ((randomb1>=0.85)&&(randomb1<0.9)){
            action3 = 17;
        }
        if ((randomb1>=0.9)&&(randomb1<0.95)){
            action3 = 18;
        }
        if ((randomb1>=0.95)){
            action3 = 19;
        }
    }
    else{
        action3 = largest(Q2,20,true_st3);
    }


    if(action3==0){
       rate3=0.05;
    }
    if(action3==1){
        rate3=0.1;
    }
    if(action3==2){
        rate3=0.15;
    }
    if(action3==3){
       rate3=0.2;
    }
    if(action3==4){
        rate3=0.25;
    }
    if(action3==5){
        rate3=0.3;
    }
    if(action3==6){
       rate3=0.35;
    }
    if(action3==7){
        rate3=0.4;
    }
    if(action3==8){
        rate3=0.45;
    }
    if(action3==9){
        rate3=0.5;
    }
    if(action3==10){
       rate3=0.55;
    }
    if(action3==11){
        rate3=0.6;
    }
    if(action3==12){
        rate3=0.65;
    }
    if(action3==13){
       rate3=0.7;
    }
    if(action3==14){
        rate3=0.75;
    }
    if(action3==15){
        rate3=0.8;
    }
    if(action3==16){
       rate3=0.85;
    }
    if(action3==17){
        rate3=0.9;
    }
    if(action3==18){
        rate3=0.95;
    }
    if(action3==19){
        rate3=1.00;
    }





    float randomc = rand()/(double)(RAND_MAX);
    float randomc1 = rand()/(double)(RAND_MAX);
    float randomc2 = rand()/(double)(RAND_MAX);


    if (randomc<epsilon){
        if (randomc1<0.05){
            action4 = 0;
        }
        if ((randomc1>=0.05)&&(randomc1<0.1)){
            action4 = 1;
        }
        if ((randomc1>=0.10)&&(randomc1<0.15)){
            action4 = 2;
        }
        if ((randomc1>=0.15)&&(randomc1<0.2)){
            action4 = 3;
        }
        if ((randomc1>=0.2)&&(randomc1<0.25)){
            action4 = 4;
        }
        if ((randomc1>=0.25)&&(randomc1<0.3)){
            action4 = 5;
        }
        if ((randomc1>=0.3)&&(randomc1<0.35)){
            action4 = 6;
        }
        if ((randomc1>=0.35)&&(randomc1<0.4)){
            action4 = 7;
        }
        if ((randomc1>=0.4)&&(randomc1<0.45)){
            action4 = 8;
        }
        if ((randomc1>=0.45)&&(randomc1<0.5)){
            action4 = 9;
        }
        if ((randomc1>=0.5)&&(randomc1<0.55)){
            action4 = 10;
        }
        if ((randomc1>=0.55)&&(randomc1<0.6)){
            action4 = 11;
        }
        if ((randomc1>=0.6)&&(randomc1<0.65)){
            action4 = 12;
        }
        if ((randomc1>=0.65)&&(randomc1<0.7)){
            action4 = 13;
        }
        if ((randomc1>=0.7)&&(randomc1<0.75)){
            action4 = 14;
        }
        if ((randomc1>=0.75)&&(randomc1<0.8)){
            action4 = 15;
        }
        if ((randomc1>=0.8)&&(randomc1<0.85)){
            action4 = 16;
        }
        if ((randomc1>=0.85)&&(randomc1<0.9)){
            action4 = 17;
        }
        if ((randomc1>=0.9)&&(randomc1<0.95)){
            action4 = 18;
        }
        if ((randomc1>=0.95)){
            action4 = 19;
        }
    }
    else{
        action4 = largest(Q3,20,true_st4);
    }


    if(action4==0){
       rate4=0.05;
    }
    if(action4==1){
        rate4=0.1;
    }
    if(action4==2){
        rate4=0.15;
    }
    if(action4==3){
       rate4=0.2;
    }
    if(action4==4){
        rate4=0.25;
    }
    if(action4==5){
        rate4=0.3;
    }
    if(action4==6){
       rate4=0.35;
    }
    if(action4==7){
        rate4=0.4;
    }
    if(action4==8){
        rate4=0.45;
    }
    if(action4==9){
        rate4=0.5;
    }
    if(action4==10){
       rate4=0.55;
    }
    if(action4==11){
        rate4=0.6;
    }
    if(action4==12){
        rate4=0.65;
    }
    if(action4==13){
       rate4=0.7;
    }
    if(action4==14){
        rate4=0.75;
    }
    if(action4==15){
        rate4=0.8;
    }
    if(action4==16){
       rate4=0.85;
    }
    if(action4==17){
        rate4=0.9;
    }
    if(action4==18){
        rate4=0.95;
    }
    if(action4==19){
        rate4=1.00;
    }


    float randomd = rand()/(double)(RAND_MAX);
    float randomd1 = rand()/(double)(RAND_MAX);
    float randomd2 = rand()/(double)(RAND_MAX);


    if (randomd<epsilon){
        if (randomd1<0.05){
            action5 = 0;
        }
        if ((randomd1>=0.05)&&(randomd1<0.1)){
            action5 = 1;
        }
        if ((randomd1>=0.10)&&(randomd1<0.15)){
            action5 = 2;
        }
        if ((randomd1>=0.15)&&(randomd1<0.2)){
            action5 = 3;
        }
        if ((randomd1>=0.2)&&(randomd1<0.25)){
            action5 = 4;
        }
        if ((randomd1>=0.25)&&(randomd1<0.3)){
            action5 = 5;
        }
        if ((randomd1>=0.3)&&(randomd1<0.35)){
            action5 = 6;
        }
        if ((randomd1>=0.35)&&(randomd1<0.4)){
            action5 = 7;
        }
        if ((randomd1>=0.4)&&(randomd1<0.45)){
            action5 = 8;
        }
        if ((randomd1>=0.45)&&(randomd1<0.5)){
            action5 = 9;
        }
        if ((randomd1>=0.5)&&(randomd1<0.55)){
            action5 = 10;
        }
        if ((randomd1>=0.55)&&(randomd1<0.6)){
            action5 = 11;
        }
        if ((randomd1>=0.6)&&(randomd1<0.65)){
            action5 = 12;
        }
        if ((randomd1>=0.65)&&(randomd1<0.7)){
            action5 = 13;
        }
        if ((randomd1>=0.7)&&(randomd1<0.75)){
            action5 = 14;
        }
        if ((randomd1>=0.75)&&(randomd1<0.8)){
            action5 = 15;
        }
        if ((randomd1>=0.8)&&(randomd1<0.85)){
            action5 = 16;
        }
        if ((randomd1>=0.85)&&(randomd1<0.9)){
            action5 = 17;
        }
        if ((randomd1>=0.9)&&(randomd1<0.95)){
            action5 = 18;
        }
        if ((randomd1>=0.95)){
            action5 = 19;
        }
    }
    else{
        action5 = largest(Q4,20,true_st5);
    }


    if(action5==0){
       rate5=0.05;
    }
    if(action5==1){
        rate5=0.1;
    }
    if(action5==2){
        rate5=0.15;
    }
    if(action5==3){
       rate5=0.2;
    }
    if(action5==4){
        rate5=0.25;
    }
    if(action5==5){
        rate5=0.3;
    }
    if(action5==6){
       rate5=0.35;
    }
    if(action5==7){
        rate5=0.4;
    }
    if(action5==8){
        rate5=0.45;
    }
    if(action5==9){
        rate5=0.5;
    }
    if(action5==10){
       rate5=0.55;
    }
    if(action5==11){
        rate5=0.6;
    }
    if(action5==12){
        rate5=0.65;
    }
    if(action5==13){
       rate5=0.7;
    }
    if(action5==14){
        rate5=0.75;
    }
    if(action5==15){
        rate5=0.8;
    }
    if(action5==16){
       rate5=0.85;
    }
    if(action5==17){
        rate5=0.9;
    }
    if(action5==18){
        rate5=0.95;
    }
    if(action5==19){
        rate5=1.00;
    }





    float randome = rand()/(double)(RAND_MAX);
    float randome1 = rand()/(double)(RAND_MAX);
    float randome2 = rand()/(double)(RAND_MAX);


    if (randome<epsilon){
        if (randome1<0.05){
            action6 = 0;
        }
        if ((randome1>=0.05)&&(randome1<0.1)){
            action6 = 1;
        }
        if ((randome1>=0.10)&&(randome1<0.15)){
            action6 = 2;
        }
        if ((randome1>=0.15)&&(randome1<0.2)){
            action6 = 3;
        }
        if ((randome1>=0.2)&&(randome1<0.25)){
            action6 = 4;
        }
        if ((randome1>=0.25)&&(randome1<0.3)){
            action6 = 5;
        }
        if ((randome1>=0.3)&&(randome1<0.35)){
            action6 = 6;
        }
        if ((randome1>=0.35)&&(randome1<0.4)){
            action6 = 7;
        }
        if ((randome1>=0.4)&&(randome1<0.45)){
            action6 = 8;
        }
        if ((randome1>=0.45)&&(randome1<0.5)){
            action6 = 9;
        }
        if ((randome1>=0.5)&&(randome1<0.55)){
            action6 = 10;
        }
        if ((randome1>=0.55)&&(randome1<0.6)){
            action6 = 11;
        }
        if ((randome1>=0.6)&&(randome1<0.65)){
            action6 = 12;
        }
        if ((randome1>=0.65)&&(randome1<0.7)){
            action6 = 13;
        }
        if ((randome1>=0.7)&&(randome1<0.75)){
            action6 = 14;
        }
        if ((randome1>=0.75)&&(randome1<0.8)){
            action6 = 15;
        }
        if ((randome1>=0.8)&&(randome1<0.85)){
            action6 = 16;
        }
        if ((randome1>=0.85)&&(randome1<0.9)){
            action6 = 17;
        }
        if ((randome1>=0.9)&&(randome1<0.95)){
            action6 = 18;
        }
        if ((randome1>=0.95)){
            action6 = 19;
        }
    }
    else{
        action6 = largest(Q5,20,true_st6);
    }


    if(action6==0){
       rate6=0.05;
    }
    if(action6==1){
        rate6=0.1;
    }
    if(action6==2){
        rate6=0.15;
    }
    if(action6==3){
       rate6=0.2;
    }
    if(action6==4){
        rate6=0.25;
    }
    if(action6==5){
        rate6=0.3;
    }
    if(action6==6){
       rate6=0.35;
    }
    if(action6==7){
        rate6=0.4;
    }
    if(action6==8){
        rate6=0.45;
    }
    if(action6==9){
        rate6=0.5;
    }
    if(action6==10){
       rate6=0.55;
    }
    if(action6==11){
        rate6=0.6;
    }
    if(action6==12){
        rate6=0.65;
    }
    if(action6==13){
       rate6=0.7;
    }
    if(action6==14){
        rate6=0.75;
    }
    if(action6==15){
        rate6=0.8;
    }
    if(action6==16){
       rate6=0.85;
    }
    if(action6==17){
        rate6=0.9;
    }
    if(action6==18){
        rate6=0.95;
    }
    if(action6==19){
        rate6=1.00;
    }

    float randomf = rand()/(double)(RAND_MAX);
    float randomf1 = rand()/(double)(RAND_MAX);
    float randomf2 = rand()/(double)(RAND_MAX);


    if (randomf<epsilon){
        if (randomf1<0.05){
            action7 = 0;
        }
        if ((randomf1>=0.05)&&(randomf1<0.1)){
            action7 = 1;
        }
        if ((randomf1>=0.10)&&(randomf1<0.15)){
            action7 = 2;
        }
        if ((randomf1>=0.15)&&(randomf1<0.2)){
            action7 = 3;
        }
        if ((randomf1>=0.2)&&(randomf1<0.25)){
            action7 = 4;
        }
        if ((randomf1>=0.25)&&(randomf1<0.3)){
            action7 = 5;
        }
        if ((randomf1>=0.3)&&(randomf1<0.35)){
            action7 = 6;
        }
        if ((randomf1>=0.35)&&(randomf1<0.4)){
            action7 = 7;
        }
        if ((randomf1>=0.4)&&(randomf1<0.45)){
            action7 = 8;
        }
        if ((randomf1>=0.45)&&(randomf1<0.5)){
            action7 = 9;
        }
        if ((randomf1>=0.5)&&(randomf1<0.55)){
            action7 = 10;
        }
        if ((randomf1>=0.55)&&(randomf1<0.6)){
            action7 = 11;
        }
        if ((randomf1>=0.6)&&(randomf1<0.65)){
            action7 = 12;
        }
        if ((randomf1>=0.65)&&(randomf1<0.7)){
            action7 = 13;
        }
        if ((randomf1>=0.7)&&(randomf1<0.75)){
            action7 = 14;
        }
        if ((randomf1>=0.75)&&(randomf1<0.8)){
            action7 = 15;
        }
        if ((randomf1>=0.8)&&(randomf1<0.85)){
            action7 = 16;
        }
        if ((randomf1>=0.85)&&(randomf1<0.9)){
            action7 = 17;
        }
        if ((randomf1>=0.9)&&(randomf1<0.95)){
            action7 = 18;
        }
        if ((randomf1>=0.95)){
            action7 = 19;
        }
    }
    else{
        action7 = largest(Q6,20,true_st7);
    }


    if(action7==0){
       rate7=0.05;
    }
    if(action7==1){
        rate7=0.1;
    }
    if(action7==2){
        rate7=0.15;
    }
    if(action7==3){
       rate7=0.2;
    }
    if(action7==4){
        rate7=0.25;
    }
    if(action7==5){
        rate7=0.3;
    }
    if(action7==6){
       rate7=0.35;
    }
    if(action7==7){
        rate7=0.4;
    }
    if(action7==8){
        rate7=0.45;
    }
    if(action7==9){
        rate7=0.5;
    }
    if(action7==10){
       rate7=0.55;
    }
    if(action7==11){
        rate7=0.6;
    }
    if(action7==12){
        rate7=0.65;
    }
    if(action7==13){
       rate7=0.7;
    }
    if(action7==14){
        rate7=0.75;
    }
    if(action7==15){
        rate7=0.8;
    }
    if(action7==16){
       rate7=0.85;
    }
    if(action7==17){
        rate7=0.9;
    }
    if(action7==18){
        rate7=0.95;
    }
    if(action7==19){
        rate7=1.00;
    }







    float randomg = rand()/(double)(RAND_MAX);
    float randomg1 = rand()/(double)(RAND_MAX);
    float randomg2 = rand()/(double)(RAND_MAX);


    if (randomg<epsilon){
        if (randomg1<0.05){
            action8 = 0;
        }
        if ((randomg1>=0.05)&&(randomg1<0.1)){
            action8 = 1;
        }
        if ((randomg1>=0.10)&&(randomg1<0.15)){
            action8 = 2;
        }
        if ((randomg1>=0.15)&&(randomg1<0.2)){
            action8 = 3;
        }
        if ((randomg1>=0.2)&&(randomg1<0.25)){
            action8 = 4;
        }
        if ((randomg1>=0.25)&&(randomg1<0.3)){
            action8 = 5;
        }
        if ((randomg1>=0.3)&&(randomg1<0.35)){
            action8 = 6;
        }
        if ((randomg1>=0.35)&&(randomg1<0.4)){
            action8 = 7;
        }
        if ((randomg1>=0.4)&&(randomg1<0.45)){
            action8 = 8;
        }
        if ((randomg1>=0.45)&&(randomg1<0.5)){
            action8 = 9;
        }
        if ((randomg1>=0.5)&&(randomg1<0.55)){
            action8 = 10;
        }
        if ((randomg1>=0.55)&&(randomg1<0.6)){
            action8 = 11;
        }
        if ((randomg1>=0.6)&&(randomg1<0.65)){
            action8 = 12;
        }
        if ((randomg1>=0.65)&&(randomg1<0.7)){
            action8 = 13;
        }
        if ((randomg1>=0.7)&&(randomg1<0.75)){
            action8 = 14;
        }
        if ((randomg1>=0.75)&&(randomg1<0.8)){
            action8 = 15;
        }
        if ((randomg1>=0.8)&&(randomg1<0.85)){
            action8 = 16;
        }
        if ((randomg1>=0.85)&&(randomg1<0.9)){
            action8 = 17;
        }
        if ((randomg1>=0.9)&&(randomg1<0.95)){
            action8 = 18;
        }
        if ((randomg1>=0.95)){
            action8 = 19;
        }
    }
    else{
        action8 = largest(Q7,20,true_st8);
    }


    if(action8==0){
       rate8=0.05;
    }
    if(action8==1){
        rate8=0.1;
    }
    if(action8==2){
        rate8=0.15;
    }
    if(action8==3){
       rate8=0.2;
    }
    if(action8==4){
        rate8=0.25;
    }
    if(action8==5){
        rate8=0.3;
    }
    if(action8==6){
       rate8=0.35;
    }
    if(action8==7){
        rate8=0.4;
    }
    if(action8==8){
        rate8=0.45;
    }
    if(action8==9){
        rate8=0.5;
    }
    if(action8==10){
       rate8=0.55;
    }
    if(action8==11){
        rate8=0.6;
    }
    if(action8==12){
        rate8=0.65;
    }
    if(action8==13){
       rate8=0.7;
    }
    if(action8==14){
        rate8=0.75;
    }
    if(action8==15){
        rate8=0.8;
    }
    if(action8==16){
       rate8=0.85;
    }
    if(action8==17){
        rate8=0.9;
    }
    if(action8==18){
        rate8=0.95;
    }
    if(action8==19){
        rate8=1.00;
    }


    float randomh = rand()/(double)(RAND_MAX);
    float randomh1 = rand()/(double)(RAND_MAX);
    float randomh2 = rand()/(double)(RAND_MAX);


    if (randomh<epsilon){
        if (randomh1<0.05){
            action9 = 0;
        }
        if ((randomh1>=0.05)&&(randomh1<0.1)){
            action9 = 1;
        }
        if ((randomh1>=0.10)&&(randomh1<0.15)){
            action9 = 2;
        }
        if ((randomh1>=0.15)&&(randomh1<0.2)){
            action9 = 3;
        }
        if ((randomh1>=0.2)&&(randomh1<0.25)){
            action9 = 4;
        }
        if ((randomh1>=0.25)&&(randomh1<0.3)){
            action9 = 5;
        }
        if ((randomh1>=0.3)&&(randomh1<0.35)){
            action9 = 6;
        }
        if ((randomh1>=0.35)&&(randomh1<0.4)){
            action9 = 7;
        }
        if ((randomh1>=0.4)&&(randomh1<0.45)){
            action9 = 8;
        }
        if ((randomh1>=0.45)&&(randomh1<0.5)){
            action9 = 9;
        }
        if ((randomh1>=0.5)&&(randomh1<0.55)){
            action9 = 10;
        }
        if ((randomh1>=0.55)&&(randomh1<0.6)){
            action9 = 11;
        }
        if ((randomh1>=0.6)&&(randomh1<0.65)){
            action9 = 12;
        }
        if ((randomh1>=0.65)&&(randomh1<0.7)){
            action9 = 13;
        }
        if ((randomh1>=0.7)&&(randomh1<0.75)){
            action9 = 14;
        }
        if ((randomh1>=0.75)&&(randomh1<0.8)){
            action9 = 15;
        }
        if ((randomh1>=0.8)&&(randomh1<0.85)){
            action9 = 16;
        }
        if ((randomh1>=0.85)&&(randomh1<0.9)){
            action9 = 17;
        }
        if ((randomh1>=0.9)&&(randomh1<0.95)){
            action9 = 18;
        }
        if ((randomh1>=0.95)){
            action9 = 19;
        }
    }
    else{
        action9 = largest(Q8,20,true_st9);
    }


    if(action9==0){
       rate9=0.05;
    }
    if(action9==1){
        rate9=0.1;
    }
    if(action9==2){
        rate9=0.15;
    }
    if(action9==3){
       rate9=0.2;
    }
    if(action9==4){
        rate9=0.25;
    }
    if(action9==5){
        rate9=0.3;
    }
    if(action9==6){
       rate9=0.35;
    }
    if(action9==7){
        rate9=0.4;
    }
    if(action9==8){
        rate9=0.45;
    }
    if(action9==9){
        rate9=0.5;
    }
    if(action9==10){
       rate9=0.55;
    }
    if(action9==11){
        rate9=0.6;
    }
    if(action9==12){
        rate9=0.65;
    }
    if(action9==13){
       rate9=0.7;
    }
    if(action9==14){
        rate9=0.75;
    }
    if(action9==15){
        rate9=0.8;
    }
    if(action9==16){
       rate9=0.85;
    }
    if(action9==17){
        rate9=0.9;
    }
    if(action9==18){
        rate9=0.95;
    }
    if(action9==19){
        rate9=1.00;
    }


    float randomi = rand()/(double)(RAND_MAX);
    float randomi1 = rand()/(double)(RAND_MAX);
    float randomi2 = rand()/(double)(RAND_MAX);


    if (randomi<epsilon){
        if (randomi1<0.05){
            action10 = 0;
        }
        if ((randomi1>=0.05)&&(randomi1<0.1)){
            action10 = 1;
        }
        if ((randomi1>=0.10)&&(randomi1<0.15)){
            action10 = 2;
        }
        if ((randomi1>=0.15)&&(randomi1<0.2)){
            action10 = 3;
        }
        if ((randomi1>=0.2)&&(randomi1<0.25)){
            action10 = 4;
        }
        if ((randomi1>=0.25)&&(randomi1<0.3)){
            action10 = 5;
        }
        if ((randomi1>=0.3)&&(randomi1<0.35)){
            action10 = 6;
        }
        if ((randomi1>=0.35)&&(randomi1<0.4)){
            action10 = 7;
        }
        if ((randomi1>=0.4)&&(randomi1<0.45)){
            action10 = 8;
        }
        if ((randomi1>=0.45)&&(randomi1<0.5)){
            action10 = 9;
        }
        if ((randomi1>=0.5)&&(randomi1<0.55)){
            action10 = 10;
        }
        if ((randomi1>=0.55)&&(randomi1<0.6)){
            action10 = 11;
        }
        if ((randomi1>=0.6)&&(randomi1<0.65)){
            action10 = 12;
        }
        if ((randomi1>=0.65)&&(randomi1<0.7)){
            action10 = 13;
        }
        if ((randomi1>=0.7)&&(randomi1<0.75)){
            action10 = 14;
        }
        if ((randomi1>=0.75)&&(randomi1<0.8)){
            action10 = 15;
        }
        if ((randomi1>=0.8)&&(randomi1<0.85)){
            action10 = 16;
        }
        if ((randomi1>=0.85)&&(randomi1<0.9)){
            action10 = 17;
        }
        if ((randomi1>=0.9)&&(randomi1<0.95)){
            action10 = 18;
        }
        if ((randomi1>=0.95)){
            action10 = 19;
        }
    }
    else{
        action10 = largest(Q9,20,true_st10);
    }


    if(action10==0){
       rate10=0.05;
    }
    if(action10==1){
        rate10=0.1;
    }
    if(action10==2){
        rate10=0.15;
    }
    if(action10==3){
       rate10=0.2;
    }
    if(action10==4){
        rate10=0.25;
    }
    if(action10==5){
        rate10=0.3;
    }
    if(action10==6){
       rate10=0.35;
    }
    if(action10==7){
        rate10=0.4;
    }
    if(action10==8){
        rate10=0.45;
    }
    if(action10==9){
        rate10=0.5;
    }
    if(action10==10){
       rate10=0.55;
    }
    if(action10==11){
        rate10=0.6;
    }
    if(action10==12){
        rate10=0.65;
    }
    if(action10==13){
       rate10=0.7;
    }
    if(action10==14){
        rate10=0.75;
    }
    if(action10==15){
        rate10=0.8;
    }
    if(action10==16){
       rate10=0.85;
    }
    if(action10==17){
        rate10=0.9;
    }
    if(action10==18){
        rate10=0.95;
    }
    if(action10==19){
        rate10=1.00;
    }


    float randomj = rand()/(double)(RAND_MAX);
    float randomj1 = rand()/(double)(RAND_MAX);
    float randomj2 = rand()/(double)(RAND_MAX);


    if (randomj<epsilon){
        if (randomj1<0.05){
            action11 = 0;
        }
        if ((randomj1>=0.05)&&(randomj1<0.1)){
            action11 = 1;
        }
        if ((randomj1>=0.10)&&(randomj1<0.15)){
            action11 = 2;
        }
        if ((randomj1>=0.15)&&(randomj1<0.2)){
            action11 = 3;
        }
        if ((randomj1>=0.2)&&(randomj1<0.25)){
            action11 = 4;
        }
        if ((randomj1>=0.25)&&(randomj1<0.3)){
            action11 = 5;
        }
        if ((randomj1>=0.3)&&(randomj1<0.35)){
            action11 = 6;
        }
        if ((randomj1>=0.35)&&(randomj1<0.4)){
            action11 = 7;
        }
        if ((randomj1>=0.4)&&(randomj1<0.45)){
            action11 = 8;
        }
        if ((randomj1>=0.45)&&(randomj1<0.5)){
            action11 = 9;
        }
        if ((randomj1>=0.5)&&(randomj1<0.55)){
            action11 = 10;
        }
        if ((randomj1>=0.55)&&(randomj1<0.6)){
            action11 = 11;
        }
        if ((randomj1>=0.6)&&(randomj1<0.65)){
            action11 = 12;
        }
        if ((randomj1>=0.65)&&(randomj1<0.7)){
            action11 = 13;
        }
        if ((randomj1>=0.7)&&(randomj1<0.75)){
            action11 = 14;
        }
        if ((randomj1>=0.75)&&(randomj1<0.8)){
            action11 = 15;
        }
        if ((randomj1>=0.8)&&(randomj1<0.85)){
            action11 = 16;
        }
        if ((randomj1>=0.85)&&(randomj1<0.9)){
            action11 = 17;
        }
        if ((randomj1>=0.9)&&(randomj1<0.95)){
            action11 = 18;
        }
        if ((randomj1>=0.95)){
            action11 = 19;
        }
    }
    else{
        action11 = largest(Q10,20,true_st11);
    }


    if(action11==0){
       rate11=0.05;
    }
    if(action11==1){
        rate11=0.1;
    }
    if(action11==2){
        rate11=0.15;
    }
    if(action11==3){
       rate11=0.2;
    }
    if(action11==4){
        rate11=0.25;
    }
    if(action11==5){
        rate11=0.3;
    }
    if(action11==6){
       rate11=0.35;
    }
    if(action11==7){
        rate11=0.4;
    }
    if(action11==8){
        rate11=0.45;
    }
    if(action11==9){
        rate11=0.5;
    }
    if(action11==10){
       rate11=0.55;
    }
    if(action11==11){
        rate11=0.6;
    }
    if(action11==12){
        rate11=0.65;
    }
    if(action11==13){
       rate11=0.7;
    }
    if(action11==14){
        rate11=0.75;
    }
    if(action11==15){
        rate11=0.8;
    }
    if(action11==16){
       rate11=0.85;
    }
    if(action11==17){
        rate11=0.9;
    }
    if(action11==18){
        rate11=0.95;
    }
    if(action11==19){
        rate11=1.00;
    }




    float randomk = rand()/(double)(RAND_MAX);
    float randomk1 = rand()/(double)(RAND_MAX);
    float randomk2 = rand()/(double)(RAND_MAX);


    if (randomk<epsilon){
        if (randomk1<0.05){
            action12 = 0;
        }
        if ((randomk1>=0.05)&&(randomk1<0.1)){
            action12 = 1;
        }
        if ((randomk1>=0.10)&&(randomk1<0.15)){
            action12 = 2;
        }
        if ((randomk1>=0.15)&&(randomk1<0.2)){
            action12 = 3;
        }
        if ((randomk1>=0.2)&&(randomk1<0.25)){
            action12 = 4;
        }
        if ((randomk1>=0.25)&&(randomk1<0.3)){
            action12 = 5;
        }
        if ((randomk1>=0.3)&&(randomk1<0.35)){
            action12 = 6;
        }
        if ((randomk1>=0.35)&&(randomk1<0.4)){
            action12 = 7;
        }
        if ((randomk1>=0.4)&&(randomk1<0.45)){
            action12 = 8;
        }
        if ((randomk1>=0.45)&&(randomk1<0.5)){
            action12 = 9;
        }
        if ((randomk1>=0.5)&&(randomk1<0.55)){
            action12 = 10;
        }
        if ((randomk1>=0.55)&&(randomk1<0.6)){
            action12 = 11;
        }
        if ((randomk1>=0.6)&&(randomk1<0.65)){
            action12 = 12;
        }
        if ((randomk1>=0.65)&&(randomk1<0.7)){
            action12 = 13;
        }
        if ((randomk1>=0.7)&&(randomk1<0.75)){
            action12 = 14;
        }
        if ((randomk1>=0.75)&&(randomk1<0.8)){
            action12 = 15;
        }
        if ((randomk1>=0.8)&&(randomk1<0.85)){
            action12 = 16;
        }
        if ((randomk1>=0.85)&&(randomk1<0.9)){
            action12 = 17;
        }
        if ((randomk1>=0.9)&&(randomk1<0.95)){
            action12 = 18;
        }
        if ((randomk1>=0.95)){
            action12 = 19;
        }
    }
    else{
        action12 = largest(Q11,20,true_st12);
    }


    if(action12==0){
       rate12=0.05;
    }
    if(action12==1){
        rate12=0.1;
    }
    if(action12==2){
        rate12=0.15;
    }
    if(action12==3){
       rate12=0.2;
    }
    if(action12==4){
        rate12=0.25;
    }
    if(action12==5){
        rate12=0.3;
    }
    if(action12==6){
       rate12=0.35;
    }
    if(action12==7){
        rate12=0.4;
    }
    if(action12==8){
        rate12=0.45;
    }
    if(action12==9){
        rate12=0.5;
    }
    if(action12==10){
       rate12=0.55;
    }
    if(action12==11){
        rate12=0.6;
    }
    if(action12==12){
        rate12=0.65;
    }
    if(action12==13){
       rate12=0.7;
    }
    if(action12==14){
        rate12=0.75;
    }
    if(action12==15){
        rate12=0.8;
    }
    if(action12==16){
       rate12=0.85;
    }
    if(action12==17){
        rate12=0.9;
    }
    if(action12==18){
        rate12=0.95;
    }
    if(action12==19){
        rate12=1.00;
    }








    pr_s1=true_s1;
    pr_s2=true_s2;
    pr_s3=true_s3;
    pr_s4=true_s4;
    pr_s5=true_s5;
    pr_s6=true_s6;
    pr_s7=true_s7;
    pr_s8=true_s8;
    pr_s9=true_s9;
    pr_s10=true_s10;
    pr_s11=true_s11;
    pr_s12=true_s12;
    pr_S=true_S;


	pr_o1=o1;
	pr_o2=o2;
	pr_o3=o3;
	pr_o4=o4;
	pr_o5=o5;
	pr_o6=o6;
	pr_o7=o7;
	pr_o8=o8;
	pr_o9=o9;
	pr_o10=o10;
	pr_o11=o11;
	pr_o12=o12;


    pr_fair1=fair1;
    pr_fair2=fair2;
    pr_fair3=fair3;
    pr_fair4=fair4;
    pr_fair5=fair5;
    pr_fair6=fair6;
    pr_fair7=fair7;
    pr_fair8=fair8;
    pr_fair9=fair9;
    pr_fair10=fair10;
    pr_fair11=fair11;
    pr_fair12=fair12;

    pr_fbk1=fbk1;
    pr_fbk2=fbk2;
    pr_fbk3=fbk3;



    for(int it=0; it<iterations; it++){
        Create_Nodes();
        //Create_Ev_Pool();

        Event_List = (Event *)NULL;
        last_inserted_event = (Event *)NULL;
        ev_access_ptr = 0;
        Channel = (NodeInfo *)NULL;
        global_time = 1.0e-16;
        Sim_Stop = FALSE;
        Initiate_Events();
        pf_total_gen_packets = 0;
        p_gen1=0;
        p_gen2=0;
        pf_total_tx_success = 0;
        sc_n1 = 0;
        sc_n2 = 0;
        sc_n3 = 0;
        sc_n4 = 0;
        sc_n5 = 0;
        sc_n6 = 0;
        sc_n7 = 0;
        sc_n8 = 0;
        sc_n9 = 0;
        sc_n10 = 0;
        sc_n11 = 0;
        sc_n12 = 0;
        coll1 = 0;
        coll2 = 0;
        coll3 = 0;
        coll4 = 0;
        coll5 = 0;
        coll6 = 0;
        coll7 = 0;
        coll8 = 0;
        coll9 = 0;
        coll10 = 0;
        coll11 = 0;
        coll12 = 0;
        coll31= 0;
        coll34= 0;
        coll42= 0;
        coll43= 0;
        new_g = 0;
        new_g1 = 0;
        new_g2 = 0;
        new_g3 = 0;
        new_g4 = 0;
        new_g5 = 0;
        new_g6 = 0;
        new_g7 = 0;
        new_g8 = 0;
        new_g9 = 0;
        new_g10 = 0;
        new_g11 = 0;
        new_g12 = 0;
        def = 0;
        def1=0;
        def2=0;
        def3=0;
        def4=0;
        drop = 0;
        drop1 =0;
        drop2=0;
        drop3=0;
        drop4=0;
        drop5=0;
        drop6=0;
        drop7=0;
        drop8=0;
        drop9=0;
        drop10=0;
        drop11=0;
        drop12=0;
        coll = 0;
        debug1 = 0;
        debug2 = 0;
        debug3 = 0;
        debug4 = 0;
        debug5 = 0;
        debug6 = 0;
        tx_flag = 0;
        n1_pkt = 0;
        n1_tx = 0;
        n1_dr=0;
        n3_pkt=0;
        n3_tx=0;
        n3_dr=0;
        n2_pkt=0;
        n2_tx=0;
        n2_dr=0;
        n4_pkt = 0;
        n4_tx = 0;
        n4_dr=0;
        n5_pkt = 0;
        n5_tx = 0;
        n5_dr=0;
        n5_pkt = 0;
        n5_tx = 0;
        n5_dr=0;
        n6_pkt = 0;
        n6_tx = 0;
        n6_dr=0;
        n7_pkt = 0;
        n7_tx = 0;
        n7_dr=0;
        n8_pkt = 0;
        n8_tx = 0;
        n8_dr=0;
        n9_pkt = 0;
        n9_tx = 0;
        n9_dr=0;
        n10_pkt = 0;
        n10_tx = 0;
        n10_dr=0;
        n11_pkt = 0;
        n11_tx = 0;
        n11_dr=0;
        n12_pkt = 0;
        n12_tx = 0;
        n12_dr=0;
        coll_ct13=0;
        coll_ct24=0;
        coll_ct34=0;
        coll_ct43=0;



        Simulate();


    //printf("s ratios: %f\t%f\tload ratios: %f\t%f\t", s_ratio1, s_ratio2, load_ratio1, load_ratio2);


    /*

            s1_vector[it]=s1_from_2;
            s2_vector[it]=(s2_from_1+s2_from_3+s2_from_4)/3;
            s3_vector[it]=(s3_from_2+s3_from_4+s3_from_5)/3;
            s4_vector[it]=0.5*(s4_from_2+s4_from_3);
            s5_vector[it]=s5_from_3;
            s6_vector[it]=(s3_from_2+s3_from_4+s3_from_5)/3;
            s7_vector[it]=0.5*(s4_from_2+s4_from_3);
            s8_vector[it]=s5_from_3;
            s1_self=s1_from_2;
            s2_self=(s2_from_1+s2_from_3+s2_from_4)/3;
            s3_self=(s3_from_2+s3_from_4+s3_from_5)/3;
            s4_self=0.5*(s4_from_2+s4_from_3);
            s5_self=s5_from_3;
            */



            s1_vector[it]=s1;
            s2_vector[it]=s2;
            s3_vector[it]=s3;
            s4_vector[it]=s4;
            s5_vector[it]=s5;
            s6_vector[it]=s6;
            s7_vector[it]=s7;
            s8_vector[it]=s8;
            s9_vector[it]=s9;
            s10_vector[it]=s10;
            s11_vector[it]=s11;
            s12_vector[it]=s12;

            //s12_vector[it]=s12;
            s21_vector[it]=s21;
            s23_vector[it]=s23;
            s24_vector[it]=s24;
            s32_vector[it]=s32;
            s34_vector[it]=s34;
            s35_vector[it]=s35;
            s42_vector[it]=s42;
            s43_vector[it]=s43;
            s53_vector[it]=s53;

    }

    /*
    for(int it=0; it<iterations; it++){
        if(count_vec[iter]==1){
            true_s1=s1_vector[it];
            true_s2=s2_vector[it];
            true_s3=s3_vector[it];
            true_s4=s4_vector[it];
            true_st1=st1_vector[it];
            true_st2=st2_vector[it];
            true_st3=st3_vector[it];
            true_st4=st4_vector[it];
        }
    }
    */

    float sum1=0.0;
    float sum2=0.0;
    float sum3=0.0;
    float sum4=0.0;
    float sum5=0.0;
    float sum6=0.0;
    float sum7=0.0;
    float sum8=0.0;
    float sum9=0.0;
    float sum10=0.0;
    float sum11=0.0;
    float sum12=0.0;
    float sum21=0.0;
    float sum23=0.0;
    float sum24=0.0;
    float sum32=0.0;
    float sum34=0.0;
    float sum35=0.0;
    float sum42=0.0;
    float sum43=0.0;
    float sum53=0.0;

    for(int it=0; it<iterations; it++){
       sum1=sum1+s1_vector[it];
       sum2=sum2+s2_vector[it];
       sum3=sum3+s3_vector[it];
       sum4=sum4+s4_vector[it];
       sum5=sum5+s5_vector[it];
       sum6=sum6+s6_vector[it];
       sum7=sum7+s7_vector[it];
       sum8=sum8+s8_vector[it];
       sum9=sum9+s9_vector[it];
       sum10=sum10+s10_vector[it];
       sum11=sum11+s11_vector[it];
       sum12=sum12+s12_vector[it];
       /*
       sum21=sum21+s21_vector[it];
       //sum12=sum12+s12_vector[it];
       sum23=sum23+s23_vector[it];
       sum24=sum24+s24_vector[it];
       sum32=sum32+s32_vector[it];
       sum34=sum34+s34_vector[it];
       sum35=sum35+s35_vector[it];
       sum42=sum42+s42_vector[it];
       sum43=sum43+s43_vector[it];
       sum53=sum53+s53_vector[it];
       */
    }

            true_s1=sum1/8.0;
            true_s2=sum2/8.0;
            true_s3=sum3/8.0;
            true_s4=sum4/8.0;
            true_s5=sum5/8.0;
            true_s6=sum6/8.0;
            true_s7=sum7/8.0;
            true_s8=sum8/8.0;
            true_s9=sum9/8.0;
            true_s10=sum10/8.0;
            true_s11=sum11/8.0;
            true_s12=sum12/8.0;

            est_s12=sum12/8.0;
            est_s21=sum21/8.0;
            est_s23=sum23/8.0;
            est_s24=sum24/8.0;
            est_s32=sum32/8.0;
            est_s34=sum34/8.0;
            est_s35=sum35/8.0;
            est_s42=sum42/8.0;
            est_s43=sum43/8.0;
            est_s53=sum53/8.0;



            true_st1=s_id;
            true_st2=s_id1;
            true_st3=s_id2;
            true_st4=s_id3;
            true_st5=s_id4;
            true_st6=s_id5;
            true_st7=s_id6;
            true_st8=s_id7;
            true_st9=s_id8;
            true_st10=s_id9;
            true_st11=s_id10;
            true_st12=s_id11;



    true_S = true_s1+true_s2+true_s3+true_s4+true_s5+true_s6+true_s7+true_s8+true_s9+true_s10+true_s11+true_s12;




    /*


    fair1= -(fabsf(true_s1-true_s2)+fabsf(true_s1-true_s5)+fabsf(true_s1-true_s7)+fabsf(true_s1-true_s3)+fabsf(true_s1-true_s8)+fabsf(true_s1-true_s12));
    fair2= -(fabsf(true_s2-true_s1)+fabsf(true_s2-true_s3)+fabsf(true_s2-true_s5)+fabsf(true_s2-true_s8)+fabsf(true_s2-true_s4)+fabsf(true_s2-true_s7)+fabsf(true_s2-true_s6)+fabsf(true_s2-true_s10)+fabsf(true_s2-true_s11));
    fair3= -(fabsf(true_s3-true_s2)+fabsf(true_s3-true_s4)+fabsf(true_s3-true_s6)+fabsf(true_s3-true_s1)+fabsf(true_s3-true_s5)+fabsf(true_s3-true_s8)+fabsf(true_s3-true_s9));
    fair4= -(fabsf(true_s4-true_s3)+fabsf(true_s4-true_s6)+fabsf(true_s4-true_s9)+fabsf(true_s4-true_s2)+fabsf(true_s4-true_s8)+fabsf(true_s4-true_s10)+fabsf(true_s4-true_s11));
    fair5= -(fabsf(true_s5-true_s1)+fabsf(true_s5-true_s2)+fabsf(true_s5-true_s7)+fabsf(true_s5-true_s3)+fabsf(true_s5-true_s8)+fabsf(true_s5-true_s12));
    fair6= -(fabsf(true_s6-true_s3)+fabsf(true_s6-true_s4)+fabsf(true_s6-true_s8)+fabsf(true_s6-true_s2)+fabsf(true_s6-true_s7)+fabsf(true_s6-true_s9)+fabsf(true_s6-true_s10)+fabsf(true_s6-true_s11));
    fair7= -(fabsf(true_s7-true_s1)+fabsf(true_s7-true_s5)+fabsf(true_s7-true_s8)+fabsf(true_s7-true_s12)+fabsf(true_s7-true_s2)+fabsf(true_s7-true_s6)+fabsf(true_s7-true_s10)+fabsf(true_s7-true_s11));
    fair8= -(fabsf(true_s8-true_s2)+fabsf(true_s8-true_s10)+fabsf(true_s8-true_s6)+fabsf(true_s8-true_s7)+fabsf(true_s8-true_s11)+fabsf(true_s8-true_s1)+fabsf(true_s8-true_s3)+fabsf(true_s8-true_s4)+fabsf(true_s8-true_s5)+fabsf(true_s8-true_s9)+fabsf(true_s8-true_s12));
    fair9= -(fabsf(true_s9-true_s4)+fabsf(true_s9-true_s10)+fabsf(true_s9-true_s12)+fabsf(true_s9-true_s3)+fabsf(true_s9-true_s6)+fabsf(true_s9-true_s7)+fabsf(true_s9-true_s8)+fabsf(true_s9-true_s11));
    fair10= -(fabsf(true_s10-true_s9)+fabsf(true_s10-true_s8)+fabsf(true_s10-true_s11)+fabsf(true_s10-true_s2)+fabsf(true_s10-true_s4)+fabsf(true_s10-true_s6)+fabsf(true_s10-true_s7)+fabsf(true_s10-true_s12));
    fair11= -(fabsf(true_s11-true_s8)+fabsf(true_s11-true_s12)+fabsf(true_s11-true_s7)+fabsf(true_s11-true_s10)+fabsf(true_s11-true_s2)+fabsf(true_s11-true_s6)+fabsf(true_s11-true_s9));
    fair12= -(fabsf(true_s12-true_s9)+fabsf(true_s12-true_s7)+fabsf(true_s12-true_s11)+fabsf(true_s12-true_s1)+fabsf(true_s12-true_s4)+fabsf(true_s12-true_s5)+fabsf(true_s12-true_s8)+fabsf(true_s12-true_s10));


    */







    /*





    fair1= -(fabsf(true_s1-true_s2)+fabsf(true_s1-true_s3)+fabsf(true_s1-true_s4)+fabsf(true_s1-true_s5)+fabsf(true_s1-true_s6)+fabsf(true_s1-true_s7)+fabsf(true_s1-true_s8)+fabsf(true_s1-true_s9)+fabsf(true_s1-true_s10)+fabsf(true_s1-true_s11)+fabsf(true_s1-true_s12));
    fair2= -(fabsf(true_s2-true_s1)+fabsf(true_s2-true_s3)+fabsf(true_s2-true_s4)+fabsf(true_s2-true_s5)+fabsf(true_s2-true_s6)+fabsf(true_s2-true_s7)+fabsf(true_s2-true_s8)+fabsf(true_s2-true_s9)+fabsf(true_s2-true_s10)+fabsf(true_s2-true_s11)+fabsf(true_s2-true_s12));
    fair3= -(fabsf(true_s3-true_s2)+fabsf(true_s3-true_s4)+fabsf(true_s3-true_s5)+fabsf(true_s3-true_s6)+fabsf(true_s3-true_s7)+fabsf(true_s3-true_s8)+fabsf(true_s3-true_s9)+fabsf(true_s3-true_s10)+fabsf(true_s3-true_s1)+fabsf(true_s3-true_s11)+fabsf(true_s3-true_s12));
    fair4= -(fabsf(true_s4-true_s1)+fabsf(true_s4-true_s2)+fabsf(true_s4-true_s3)+fabsf(true_s4-true_s5)+fabsf(true_s4-true_s6)+fabsf(true_s4-true_s7)+fabsf(true_s4-true_s8)+fabsf(true_s4-true_s9)+fabsf(true_s4-true_s10)+fabsf(true_s4-true_s11)+fabsf(true_s4-true_s12));
    fair5= -(fabsf(true_s5-true_s1)+fabsf(true_s5-true_s2)+fabsf(true_s5-true_s3)+fabsf(true_s5-true_s4)+fabsf(true_s5-true_s6)+fabsf(true_s5-true_s7)+fabsf(true_s5-true_s8)+fabsf(true_s5-true_s9)+fabsf(true_s5-true_s10)+fabsf(true_s5-true_s11)+fabsf(true_s5-true_s12));
    fair6= -(fabsf(true_s6-true_s1)+fabsf(true_s6-true_s2)+fabsf(true_s6-true_s3)+fabsf(true_s6-true_s4)+fabsf(true_s5-true_s6)+fabsf(true_s6-true_s7)+fabsf(true_s6-true_s8)+fabsf(true_s6-true_s9)+fabsf(true_s6-true_s10)+fabsf(true_s6-true_s11)+fabsf(true_s6-true_s12));
    fair7= -(fabsf(true_s7-true_s1)+fabsf(true_s7-true_s2)+fabsf(true_s7-true_s3)+fabsf(true_s7-true_s4)+fabsf(true_s7-true_s5)+fabsf(true_s7-true_s6)+fabsf(true_s7-true_s8)+fabsf(true_s7-true_s9)+fabsf(true_s7-true_s10)+fabsf(true_s7-true_s11)+fabsf(true_s7-true_s12));
    fair8= -(fabsf(true_s8-true_s1)+fabsf(true_s8-true_s2)+fabsf(true_s8-true_s3)+fabsf(true_s8-true_s4)+fabsf(true_s8-true_s5)+fabsf(true_s8-true_s6)+fabsf(true_s8-true_s7)+fabsf(true_s8-true_s9)+fabsf(true_s8-true_s10)+fabsf(true_s8-true_s11)+fabsf(true_s8-true_s12));
    fair9= -(fabsf(true_s9-true_s1)+fabsf(true_s9-true_s2)+fabsf(true_s9-true_s3)+fabsf(true_s9-true_s4)+fabsf(true_s9-true_s5)+fabsf(true_s9-true_s6)+fabsf(true_s9-true_s7)+fabsf(true_s9-true_s8)+fabsf(true_s9-true_s10)+fabsf(true_s9-true_s11)+fabsf(true_s9-true_s12));
    fair10= -(fabsf(true_s10-true_s1)+fabsf(true_s10-true_s2)+fabsf(true_s10-true_s3)+fabsf(true_s10-true_s4)+fabsf(true_s10-true_s5)+fabsf(true_s10-true_s6)+fabsf(true_s10-true_s7)+fabsf(true_s10-true_s8)+fabsf(true_s10-true_s9)+fabsf(true_s10-true_s11)+fabsf(true_s10-true_s12));
    fair11= -(fabsf(true_s11-true_s1)+fabsf(true_s11-true_s2)+fabsf(true_s11-true_s3)+fabsf(true_s11-true_s4)+fabsf(true_s11-true_s5)+fabsf(true_s11-true_s6)+fabsf(true_s11-true_s7)+fabsf(true_s11-true_s8)+fabsf(true_s11-true_s9)+fabsf(true_s11-true_s10)+fabsf(true_s11-true_s12));
    fair12= -(fabsf(true_s12-true_s1)+fabsf(true_s12-true_s2)+fabsf(true_s12-true_s3)+fabsf(true_s12-true_s4)+fabsf(true_s12-true_s5)+fabsf(true_s12-true_s6)+fabsf(true_s12-true_s7)+fabsf(true_s12-true_s8)+fabsf(true_s12-true_s9)+fabsf(true_s12-true_s11)+fabsf(true_s12-true_s10));





    */



    fair1= -(fabsf(true_s1-true_s2)+fabsf(true_s1-true_s5)+fabsf(true_s1-true_s7));
    fair2= -(fabsf(true_s2-true_s1)+fabsf(true_s2-true_s3)+fabsf(true_s2-true_s5)+fabsf(true_s2-true_s8));
    fair3= -(fabsf(true_s3-true_s2)+fabsf(true_s3-true_s4)+fabsf(true_s3-true_s6));
    fair4= -(fabsf(true_s4-true_s3)+fabsf(true_s4-true_s6)+fabsf(true_s4-true_s9));
    fair5= -(fabsf(true_s5-true_s1)+fabsf(true_s5-true_s2)+fabsf(true_s5-true_s7));
    fair6= -(fabsf(true_s6-true_s3)+fabsf(true_s6-true_s4)+fabsf(true_s6-true_s8));
    fair7= -(fabsf(true_s7-true_s1)+fabsf(true_s7-true_s5)+fabsf(true_s7-true_s8)+fabsf(true_s7-true_s12));
    fair8= -(fabsf(true_s8-true_s2)+fabsf(true_s8-true_s6)+fabsf(true_s8-true_s7)+fabsf(true_s8-true_s10)+fabsf(true_s8-true_s11));
    fair9= -(fabsf(true_s9-true_s8)+fabsf(true_s9-true_s10)+fabsf(true_s9-true_s11));
    fair10= -(fabsf(true_s10-true_s8)+fabsf(true_s10-true_s9)+fabsf(true_s10-true_s11));
    fair11= -(fabsf(true_s11-true_s7)+fabsf(true_s11-true_s8)+fabsf(true_s11-true_s10)+fabsf(true_s11-true_s12));
    fair12= -(fabsf(true_s12-true_s7)+fabsf(true_s12-true_s9)+fabsf(true_s12-true_s11));







    /*



    fair1= -(fabsf(true_s1-true_s2)+fabsf(true_s1-true_s3)+fabsf(true_s1-true_s4));
    fair2= -(fabsf(true_s2-true_s1)+fabsf(true_s2-true_s3)+fabsf(true_s2-true_s4));
    fair3= -(fabsf(true_s3-true_s2)+fabsf(true_s3-true_s4)+fabsf(true_s3-true_s5));
    fair4= -(fabsf(true_s4-true_s3)+fabsf(true_s4-true_s5)+fabsf(true_s4-true_s6));
    fair5= -(fabsf(true_s5-true_s4)+fabsf(true_s5-true_s6)+fabsf(true_s5-true_s7));
    fair6= -(fabsf(true_s6-true_s5)+fabsf(true_s6-true_s7)+fabsf(true_s6-true_s8));
    fair7= -(fabsf(true_s7-true_s6)+fabsf(true_s7-true_s8)+fabsf(true_s7-true_s9));
    fair8= -(fabsf(true_s8-true_s7)+fabsf(true_s8-true_s9)+fabsf(true_s8-true_s10));
    fair9= -(fabsf(true_s9-true_s8)+fabsf(true_s9-true_s10)+fabsf(true_s9-true_s11));
    fair10= -(fabsf(true_s10-true_s9)+fabsf(true_s10-true_s11)+fabsf(true_s10-true_s12));
    fair11= -(fabsf(true_s11-true_s1)+fabsf(true_s11-true_s10)+fabsf(true_s11-true_s12));
    fair12= -(fabsf(true_s12-true_s1)+fabsf(true_s12-true_s2)+fabsf(true_s12-true_s11));

    */




    fair31 = -fabsf(true_s1-true_s3);
    fair32 = -fabsf(true_s2-true_s3);


    fbk1=(true_s3);
    fbk2=(true_s3);
    fbk3=0.5*(true_s1+true_s2);


    /*
    fbk1=0.5*(true_s1+true_s3);
    fbk2=0.5*(true_s2+true_s3);
    fbk3=0.33*(true_s1+true_s2+true_s3);
    */



    /*

    //One hop info
	o1 = true_s1+true_s2+true_s5+true_s7;
	o2 = true_s1+true_s2+true_s3+true_s5+true_s8;
	o3 = true_s2+true_s3+true_s4+true_s6;
	o4 = true_s3+true_s4+true_s6+true_s9;
	o5 = true_s1+true_s2+true_s5+true_s7;
	o6 = true_s3+true_s4+true_s6+true_s8;
	o7 = true_s1+true_s12+true_s5+true_s7+true_s8;
	o8 = true_s2+true_s6+true_s7+true_s8+true_s10+true_s11;
	o9 = true_s4+true_s9+true_s10+true_s12;
	o10 = true_s9+true_s11+true_s8+true_s10;
	o11 = true_s10+true_s12+true_s7+true_s11+true_s8;
	o12 = true_s7+true_s11+true_s12+true_s9;



	*/






	//Two-hop info

	o1 = true_s1+true_s2+true_s5+true_s7+true_s3+true_s12+true_s8;
	o2 = true_s1+true_s2+true_s3+true_s5+true_s8+true_s7+true_s4+true_s6+true_s10+true_s11;
	o3 = true_s2+true_s3+true_s4+true_s6+true_s5+true_s8+true_s1+true_s9;
	o4 = true_s3+true_s4+true_s6+true_s9+true_s8+true_s2+true_s10+true_s12;
	o5 = true_s1+true_s2+true_s5+true_s7+true_s3+true_s8+true_s12;
	o6 = true_s3+true_s4+true_s6+true_s8+true_s2+true_s9+true_s7+true_s10+true_s11;
	o7 = true_s1+true_s12+true_s5+true_s7+true_s8+true_s2+true_s6+true_s10+true_s11;
	o8 = true_S;
	o9 = true_s4+true_s9+true_s10+true_s12+true_s3+true_s6+true_s11+true_s8+true_s7;
	o10 = true_s9+true_s11+true_s8+true_s10+true_s2+true_s6+true_s7+true_s4+true_s12;
	o11 = true_s10+true_s12+true_s7+true_s11+true_s8+true_s2+true_s6+true_s9;
	o12 = true_s7+true_s11+true_s12+true_s9+true_s1+true_s8+true_s5+true_s10+true_s4;












    if(((o1-pr_o1)>=0.005)&&(fair1-pr_fair1)>=-0.01){
        rn1=50.0;
    }if(((o1-pr_o1)>=0.005)&&(fair1-pr_fair1)<-0.01){
        rn1=-50.0;
    }if(((o1-pr_o1)<0.005)&&(fair1-pr_fair1)>=-0.01){
        rn1=-50.0;
    }if(((o1-pr_o1)<0.005)&&(fair1-pr_fair1)<-0.01){
        rn1=-50.0;
    }
    if(((o2-pr_o2)>=0.005)&&(fair2-pr_fair2)>=-0.01){
        rn2=50.0;
    }if(((o2-pr_o2)>=0.005)&&(fair2-pr_fair2)<-0.01){
        rn2=-50.0;
    }if(((o2-pr_o2)<0.005)&&(fair2-pr_fair2)>=-0.01){
        rn2=-50.0;
    }if(((o2-pr_o2)<0.005)&&(fair2-pr_fair2)<-0.01){
        rn2=-50.0;
    }
    if(((o3-pr_o3)>=0.005)&&(fair3-pr_fair3)>=-0.01){
        rn3=50.0;
    }if(((o3-pr_o3)>=0.005)&&(fair3-pr_fair3)<-0.01){
        rn3=-50.0;
    }if(((o3-pr_o3)<0.005)&&(fair3-pr_fair3)>=-0.01){
        rn3=-50.0;
    }if(((o3-pr_o3)<0.005)&&(fair3-pr_fair3)<-0.01){
        rn3=-50.0;
    }
    if(((o4-pr_o4)>=0.005)&&(fair4-pr_fair4)>=-0.01){
        rn4=50.0;
    }if(((o4-pr_o4)>=0.005)&&(fair4-pr_fair4)<-0.01){
        rn4=-50.0;
    }if(((o4-pr_o4)<0.005)&&(fair4-pr_fair4)>=-0.01){
        rn4=-50.0;
    }if(((o4-pr_o4)<0.005)&&(fair4-pr_fair4)<-0.01){
        rn4=-50.0;
    }
    if(((o5-pr_o5)>=0.005)&&(fair5-pr_fair5)>=-0.01){
        rn5=50.0;
    }if(((o5-pr_o5)>=0.005)&&(fair5-pr_fair5)<-0.01){
        rn5=-50.0;
    }if(((o5-pr_o5)<0.005)&&(fair5-pr_fair5)>=-0.01){
        rn5=-50.0;
    }if(((o5-pr_o5)<0.005)&&(fair5-pr_fair5)<-0.01){
        rn5=-50.0;
    }
    if(((o6-pr_o6)>=0.005)&&(fair6-pr_fair6)>=-0.01){
        rn6=50.0;
    }if(((o6-pr_o6)>=0.005)&&(fair6-pr_fair6)<-0.01){
        rn6=-50.0;
    }if(((o6-pr_o6)<0.005)&&(fair6-pr_fair6)>=-0.01){
        rn6=-50.0;
    }if(((o6-pr_o6)<0.005)&&(fair6-pr_fair6)<-0.01){
        rn6=-50.0;
    }
    if(((o7-pr_o7)>=0.005)&&(fair7-pr_fair7)>=-0.01){
        rn7=50.0;
    }if(((o7-pr_o7)>=0.005)&&(fair7-pr_fair7)<-0.01){
        rn7=-50.0;
    }if(((o7-pr_o7)<0.005)&&(fair7-pr_fair7)>=-0.01){
        rn7=-50.0;
    }if(((o7-pr_o7)<0.005)&&(fair7-pr_fair7)<-0.01){
        rn7=-50.0;
    }
    if(((o8-pr_o8)>=0.005)&&(fair8-pr_fair8)>=-0.01){
        rn8=50.0;
    }if(((o8-pr_o8)>=0.005)&&(fair8-pr_fair8)<-0.01){
        rn8=-50.0;
    }if(((o8-pr_o8)<0.005)&&(fair8-pr_fair8)>=-0.01){
        rn8=-50.0;
    }if(((o8-pr_o8)<0.005)&&(fair8-pr_fair8)<-0.01){
        rn8=-50.0;
    }
    if(((o9-pr_o9)>=0.005)&&(fair9-pr_fair9)>=-0.01){
        rn9=50.0;
    }if(((o9-pr_o9)>=0.005)&&(fair9-pr_fair9)<-0.01){
        rn9=-50.0;
    }if(((o9-pr_o9)<0.005)&&(fair9-pr_fair9)>=-0.01){
        rn9=-50.0;
    }if(((o9-pr_o9)<0.005)&&(fair9-pr_fair9)<-0.01){
        rn9=-50.0;
    }
    if(((o10-pr_o10)>=0.005)&&(fair10-pr_fair10)>=-0.01){
        rn10=50.0;
    }if(((o10-pr_o10)>=0.005)&&(fair10-pr_fair10)<-0.01){
        rn10=-50.0;
    }if(((o10-pr_o10)<0.005)&&(fair10-pr_fair10)>=-0.01){
        rn10=-50.0;
    }if(((o10-pr_o10)<0.005)&&(fair10-pr_fair10)<-0.01){
        rn10=-50.0;
    }
    if(((o11-pr_o11)>=0.005)&&(fair11-pr_fair11)>=-0.01){
        rn11=50.0;
    }if(((o11-pr_o11)>=0.005)&&(fair11-pr_fair11)<-0.01){
        rn11=-50.0;
    }if(((o11-pr_o11)<0.005)&&(fair11-pr_fair11)>=-0.01){
        rn11=-50.0;
    }if(((o11-pr_o11)<0.005)&&(fair11-pr_fair11)<-0.01){
        rn11=-50.0;
    }
    if(((o12-pr_o12)>=0.005)&&(fair12-pr_fair12)>=-0.01){
        rn12=50.0;
    }if(((o12-pr_o12)>=0.005)&&(fair12-pr_fair12)<-0.01){
        rn12=-50.0;
    }if(((o12-pr_o12)<0.005)&&(fair12-pr_fair12)>=-0.01){
        rn12=-50.0;
    }if(((o12-pr_o12)<0.005)&&(fair12-pr_fair12)<-0.01){
        rn12=-50.0;
    }








    /*


    if(((true_S-pr_S)>=0.005)&&(fair1-pr_fair1)>=-0.005){
        rn1=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair1-pr_fair1)<-0.005){
        rn1=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair1-pr_fair1)>=-0.005){
        rn1=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair1-pr_fair1)<-0.005){
        rn1=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair2-pr_fair2)>=-0.005){
        rn2=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair2-pr_fair2)<-0.005){
        rn2=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair2-pr_fair2)>=-0.005){
        rn2=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair2-pr_fair2)<-0.005){
        rn2=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair3-pr_fair3)>=-0.005){
        rn3=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair3-pr_fair3)<-0.005){
        rn3=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair3-pr_fair3)>=-0.005){
        rn3=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair3-pr_fair3)<-0.005){
        rn3=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair4-pr_fair4)>=-0.005){
        rn4=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair4-pr_fair4)<-0.005){
        rn4=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair4-pr_fair4)>=-0.005){
        rn4=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair4-pr_fair4)<-0.005){
        rn4=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair5-pr_fair5)>=-0.005){
        rn5=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair5-pr_fair5)<-0.005){
        rn5=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair5-pr_fair5)>=-0.005){
        rn5=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair5-pr_fair5)<-0.005){
        rn5=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair6-pr_fair6)>=-0.005){
        rn6=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair6-pr_fair6)<-0.005){
        rn6=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair6-pr_fair6)>=-0.005){
        rn6=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair6-pr_fair6)<-0.005){
        rn6=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair7-pr_fair7)>=-0.005){
        rn7=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair7-pr_fair7)<-0.005){
        rn7=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair7-pr_fair7)>=-0.005){
        rn7=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair7-pr_fair7)<-0.005){
        rn7=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair8-pr_fair8)>=-0.005){
        rn8=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair8-pr_fair8)<-0.005){
        rn8=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair8-pr_fair8)>=-0.005){
        rn8=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair8-pr_fair8)<-0.005){
        rn8=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair9-pr_fair9)>=-0.005){
        rn9=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair9-pr_fair9)<-0.005){
        rn9=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair9-pr_fair9)>=-0.005){
        rn9=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair9-pr_fair9)<-0.005){
        rn9=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair10-pr_fair10)>=-0.005){
        rn10=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair10-pr_fair10)<-0.005){
        rn10=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair10-pr_fair10)>=-0.005){
        rn10=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair10-pr_fair10)<-0.005){
        rn10=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair11-pr_fair11)>=-0.005){
        rn11=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair11-pr_fair11)<-0.005){
        rn11=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair11-pr_fair11)>=-0.005){
        rn11=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair11-pr_fair11)<-0.005){
        rn11=-50.0;
    }
    if(((true_S-pr_S)>=0.005)&&(fair12-pr_fair12)>=-0.005){
        rn12=50.0;
    }if(((true_S-pr_S)>=0.005)&&(fair12-pr_fair12)<-0.005){
        rn12=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair12-pr_fair12)>=-0.005){
        rn12=-50.0;
    }if(((true_S-pr_S)<0.005)&&(fair12-pr_fair12)<-0.005){
        rn12=-50.0;
    }

    */








    if(true_s1==0.0){
        rn1=-100.0;
    }
    if(true_s2==0.0){
        rn2=-100.0;
    }
    if(true_s3==0.0){
        rn3=-100.0;
    }
    if(true_s4==0.0){
        rn4=-100.0;
    }
    if(true_s5==0.0){
        rn5=-100.0;
    }
    if(true_s6==0.0){
        rn6=-100.0;
    }
    if(true_s7==0.0){
        rn7=-100.0;
    }
    if(true_s8==0.0){
        rn8=-100.0;
    }
    if(true_s9==0.0){
        rn9=-100.0;
    }
    if(true_s10==0.0){
        rn10=-100.0;
    }
    if(true_s11==0.0){
        rn11=-100.0;
    }
    if(true_s12==0.0){
        rn12=-100.0;
    }







        delta1 = ((float)rn1+(float)gamma*(float)maximum(Q,20,true_st1))-(float)Q[prev_s_id][action1];
        delta2 = ((float)rn2+(float)gamma*(float)maximum(Q1,20,true_st2))-(float)Q1[prev_s_id1][action2];
        delta3 = ((float)rn3+(float)gamma*(float)maximum(Q2,20,true_st3))-(float)Q2[prev_s_id2][action3];
        delta4 = ((float)rn4+(float)gamma*(float)maximum(Q3,20,true_st4))-(float)Q3[prev_s_id3][action4];
        delta5 = ((float)rn5+(float)gamma*(float)maximum(Q4,20,true_st5))-(float)Q4[prev_s_id4][action5];
        delta6 = ((float)rn6+(float)gamma*(float)maximum(Q5,20,true_st6))-(float)Q5[prev_s_id5][action6];
        delta7 = ((float)rn7+(float)gamma*(float)maximum(Q6,20,true_st7))-(float)Q6[prev_s_id6][action7];
        delta8 = ((float)rn8+(float)gamma*(float)maximum(Q7,20,true_st8))-(float)Q7[prev_s_id7][action8];
        delta9 = ((float)rn9+(float)gamma*(float)maximum(Q8,20,true_st9))-(float)Q8[prev_s_id8][action9];
        delta10 = ((float)rn10+(float)gamma*(float)maximum(Q9,20,true_st10))-(float)Q9[prev_s_id9][action10];
        delta11 = ((float)rn11+(float)gamma*(float)maximum(Q10,20,true_st11))-(float)Q10[prev_s_id10][action11];
        delta12 = ((float)rn12+(float)gamma*(float)maximum(Q11,20,true_st12))-(float)Q11[prev_s_id11][action12];

        d1 = ((float)rn1+(float)gamma*(float)maximum(Q,19,true_st1));
        d2 = ((float)rn2+(float)gamma*(float)maximum(Q1,19,true_st2));
        d3 = ((float)rn3+(float)gamma*(float)maximum(Q2,19,true_st3));


        if(delta1>0.0){
            Q[prev_s_id][action1]=Q[prev_s_id][action1]+alpha*delta1;
        }else{
            Q[prev_s_id][action1]=Q[prev_s_id][action1]+beta*delta1;
        }


        if(delta2>0.0){
            Q1[prev_s_id1][action2]=Q1[prev_s_id1][action2]+alpha*delta2;
        }else{
            Q1[prev_s_id1][action2]=Q1[prev_s_id1][action2]+beta*delta2;
        }


        if(delta3>0.0){
            Q2[prev_s_id2][action3]=Q2[prev_s_id2][action3]+alpha*delta3;
        }else{
            Q2[prev_s_id2][action3]=Q2[prev_s_id2][action3]+beta*delta3;
        }

        if(delta4>0.0){
            Q3[prev_s_id3][action4]=Q3[prev_s_id3][action4]+alpha*delta4;
        }else{
            Q3[prev_s_id3][action4]=Q3[prev_s_id3][action4]+beta*delta4;
        }

        if(delta5>0.0){
            Q4[prev_s_id4][action5]=Q4[prev_s_id4][action5]+alpha*delta5;
        }else{
            Q4[prev_s_id4][action5]=Q4[prev_s_id4][action5]+beta*delta5;
        }
        if(delta6>0.0){
            Q5[prev_s_id5][action6]=Q5[prev_s_id5][action6]+alpha*delta6;
        }else{
            Q5[prev_s_id5][action6]=Q5[prev_s_id5][action6]+beta*delta6;
        }
        if(delta7>0.0){
            Q6[prev_s_id6][action7]=Q6[prev_s_id6][action7]+alpha*delta7;
        }else{
            Q6[prev_s_id6][action7]=Q6[prev_s_id6][action7]+beta*delta7;
        }
        if(delta8>0.0){
            Q7[prev_s_id7][action8]=Q7[prev_s_id7][action8]+alpha*delta8;
        }else{
            Q7[prev_s_id7][action8]=Q7[prev_s_id7][action8]+beta*delta8;
        }
        if(delta9>0.0){
            Q8[prev_s_id8][action9]=Q8[prev_s_id8][action9]+alpha*delta9;
        }else{
            Q8[prev_s_id8][action9]=Q8[prev_s_id8][action9]+beta*delta9;
        }
        if(delta10>0.0){
            Q9[prev_s_id9][action10]=Q9[prev_s_id9][action10]+alpha*delta10;
        }else{
            Q9[prev_s_id9][action10]=Q9[prev_s_id9][action10]+beta*delta10;
        }
        if(delta11>0.0){
            Q10[prev_s_id10][action11]=Q10[prev_s_id10][action11]+alpha*delta11;
        }else{
            Q10[prev_s_id10][action11]=Q10[prev_s_id10][action11]+beta*delta11;
        }
        if(delta12>0.0){
            Q11[prev_s_id11][action12]=Q11[prev_s_id11][action12]+alpha*delta12;
        }else{
            Q11[prev_s_id11][action12]=Q11[prev_s_id11][action12]+beta*delta12;
        }



        //Q[prev_s_id][action]=(1-(float)alpha)*(float)Q[prev_s_id][action]+(float)alpha*((float)reward+(float)gamma*(float)maximum(Q,3,s_id));


        //fprintf(fpt,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",act_g1,act_g2, act_g3, act_g4, true_s1, true_s2, true_s3, true_s4, true_S);

        fprintf(fpt,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",true_s1, true_s2, true_s3, true_s4,true_s5, true_s6, true_s7, true_s8,true_s9, true_s10, true_s11, true_s12, fair1, fair2, fair3, fair4, fair5, fair6, fair7, fair8, fair9, fair10, fair11, fair12, rn1);
        fprintf(fpt1,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",true_st1, true_st2, true_st3, true_st4,true_st5, true_st6, true_st7, true_st8,true_st9, true_st10, true_st11, true_st12, action1, action2, action3, action4, action5, action6, action7, action8, action9, action10, action11, action12, rn1, rn2, rn3, rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn11,rn12);




    printf("\n \nSimulation Ends\n\n");
    printf("Epoch ID: %d\n", epoch_id);
    printf("Total Success : %d\n",pf_total_tx_success);
    printf("Total Generated Packets : %d\n",pf_total_gen_packets);
    printf("Total Dropped Packets for busy channel: %d\n",def);
    printf("Total Dropped Packets (intentional): %d\n",drop);
    printf("Total Collided Packets : %d\n",coll);
    printf("Node 1 load: %f\t Node 2 load: %f\n",small_g1,small_g2);
    //printf("Network Load : %f\n",Capital_G);
    printf("Actual Network Load : %f\n",act_g);
    printf("Actual Node 1 Load : %f\n",act_g1);
    printf("Actual Node 2 Load : %f\n",act_g2);
    printf("Actual Node 3 Load : %f\n",act_g3);
    printf("Actual Node 4 Load : %f\n",act_g4);
    printf("Global Time : %f\n",(float)global_time);
    printf("Simulation Time: %d\n",Simulation_Duration);
    printf("Packet Duration : %d\n",packet_duration);
    //printf("Node1 events: %d\t Node2 events: %d\t ",debug1,debug2);
    printf("Node1 Transmit actions:%d\t Node1 Drop actions:%d\nNode2 Transmit actions: %d,\t Node 2 Drop Actions: %d\n",n1_tx,n1_dr,n2_tx,n2_dr);
    printf("Success ratio : %f\n",(act_g)*(double)pf_total_tx_success/(double)(new_g+0.000001));
    printf("Resulting Capital_S = %10.10f\n", (true_S));
    printf("s1:%f,  s2:%f,  s3:%f,   s4:%f,  s5:%f,  s6:%f\n   s7:%f,  s8:%f,  s9:%f,  s10:%f,   s11:%f,  s12:%f\n",true_s1,true_s2, true_s3, true_s4, true_s5, true_s6, true_s7, true_s8, true_s9, true_s10, true_s11, true_s12);
    //printf("s21:%f,  s23:%f\n",s21,s23);
    printf("a1:%d,  a2:%d,  a3:%d, a4:%d,  a5:%d\n",action1,action2, action3, action4, action5);

    //printf("%d\n%d\n%f\n%d\n%d\n",debug5, debug6, current_throughput1, pf_total_tx_success, new_g);
    //printf("Success:%f\tFailure:%f\t%f\tInter-collision:%f\n", current_throughput, pf1,pf2, prob_inter);
    //printf("Node1 success: %d\t Node2 success: %d\n",sc_n1,sc_n2);
    printf("Rate1:%f\tRate2:%f\tRate3:%f\tRate4:%f\tRate5:%f\tRate6:%f\nRate7:%f\tRate8:%f\tRate9:%f\tRate10:%f\tRate11:%f\tRate12:%f\n",rate1,rate2,rate3, rate4, rate5, rate6, rate7, rate8, rate9, rate10, rate11, rate12);
    //printf("%f,%f", alpha, beta);

    }


    printf("%d\t%d\n",cnt1,cnt2);

    //printf("%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\t%f\n%f\n",Q[10][0], Q[10][1],Q[11][0], Q[11][1], Q[12][0], Q[12][1],Q[13][0], Q[13][1],Q[14][0], Q[14][1],Q[15][0], Q[15][1], Q[16][0], Q[16][1],Q[17][0],Q[17][1],Q[18][0], Q[18][1],Q[19][0],Q[19][1], reward);
    fclose(fpt);
    fclose(fpt1);

}
