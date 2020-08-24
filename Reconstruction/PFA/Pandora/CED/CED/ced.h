/* "C" event display.
 * Main part. 
 *
 * Alexey Zhelezov, DESY/ITEP, 2005 */

/*
 * This file is internal. It must not be
 * included into enduser application.
 */

#ifndef __CED_H
#define __CED_H

#include "ced_cli.h"

//#ifdef __cplusplus
// extern "C" {
//#endif
		

//char trusted_hosts[50];
//extern static char testchar;

typedef void (*ced_draw_cb)(void *data);


/*
 * Register new element type. Order is important!
 * Appropriate code must be defined in both
 * client and server parts.
 *
 *  item_size - size of one item in bytes.
 *  draw_func - function to call to draw one item,
 *              called from ced_do_draw_event()
 *              not used on client side.
 */
unsigned ced_register_element(unsigned item_size,ced_draw_cb draw_func);

/*
 * To be called from element functions
 * on client side.
 * Return allocated space for one item
 * with size, specified by ced_register_element()
 *
 * Example: assume struct Dummy { int i; }; is item.
 *
 *          DummyID=ced_register_element(sizeof(struct Dummy),0);
 *          ...
 *          struct Dummy *item=(struct Dummy *)ced_add_element(DummyID);
 *            item->i=0;
 * This will add one Dummy item to the event
 */
void *ced_add(unsigned id);

/*
 * To be called in paint function
 *
 * It calls user defined functions for
 * each item of all elements types.
 */
void ced_do_draw_event(void);

/*
 * Server side function.
 * Must be used to process all incoming
 * messages from client.
 *
 * It return positive value when
 * new event must be drawn.
 *
 * Example:
 *      glut_tcp_server(7285,my_process_input)
 *
 *      my_process_input(x){
 *        if(ced_process_input(x)>0)
 *          <do redraw>
 */
int ced_process_input(void *data);

//------------
void addLayerDescriptionToMenu(int,char *);//glced.c
void selectFromMenu(int id);//glced.c
void toggleHelpWindow(void);//glced.c
void updateLayerEntryInPopupMenu(int); //glced.c
int buildMenuPopup(void);//glced.c


#define VERSION_CONFIG 3
struct CEDsettings{
    bool trans;         //grid or surface view
    bool persp;         //perspectivic view or flat projection 
    bool antia;         //anti aliasing
    bool light;         //light source 
    bool picking_highlight; //marker at picking position
    double detector_trans[NUMBER_DETECTOR_LAYER];
    double detector_cut_angle[NUMBER_DETECTOR_LAYER];
    double detector_cut_z[NUMBER_DETECTOR_LAYER];
    bool detector_picking;
//    double cut_angle; //deprecated
//    double trans_value; //deprecated
//    double z_cutting; //deprecated
    bool layer[CED_MAX_LAYER];
    bool phi_projection;
    bool z_projection;
    double view[3];
    double va; //vertical angle of view
    double ha; //horionzional angle of view
    bool fixed_view;
    int win_h; //height of the window (pixel)
    int win_w; //wight of the window (pixel)
    double zoom;
    double fisheye_alpha;
    double world_size;
    double fisheye_world_size;
    double bgcolor[4];
    bool show_axes;
    bool fps;
    double screenshot_sections;
    int font; //size of text (menu, shortcuts, text in ced window)
    bool autoshot; // If true, generate screencapture in every new event
    int autoshot_scale; // If true, generate screencapture in every new event
};

/*
//important: 
//          - sum of all layers must be smaler than max_layer!
//          - number_popup_layer must be smaler than number_data_layer

#define CED_MAX_LAYER       100 
#define NUMBER_POPUP_LAYER      20
#define NUMBER_DATA_LAYER       25
#define NUMBER_DETECTOR_LAYER   20


#define CED_MAX_LAYER_CHAR 400
*/

//-------------

//#ifdef __cplusplus
// }
//#endif

	

#endif /* __CED_H  */
