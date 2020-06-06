/* "C" event display.
 * Enduser accessable API.
 *
 * Alexey Zhelezov, DESY/ITEP, 2005 
 */

#ifndef __CED_CLI_H
#define __CED_CLI_H

#include <ced_config.h>


//important:
//          - sum of all layers must be smaler than max_layer!
//          - number_popup_layer must be smaler than number_data_layer

//#define CED_MAX_LAYER       100
//#define NUMBER_POPUP_LAYER      20
//#define NUMBER_DATA_LAYER       25
//#define NUMBER_DETECTOR_LAYER   20
//
//
//#define CED_MAX_LAYER_CHAR 400


#ifdef __cplusplus
extern "C" {
#endif

/*
 * This is the first function to call (before any other).
 *
 *  host - host with CED (must be "localhost")
 *  port - server port number (let say 7285 :)
 *
 * NOTE: ced_register_elements() must be called
 *       separately.
 */
void ced_client_init(const char *host,unsigned short port);

/*
 * Cancel current event output. So, all elements
 * queued will be discarded.
 *
 * Good to call at the begining of every event processing.
 */
void ced_new_event(void);


/*
 * This function really attempt to display event in CED.
 * When CED is not available, this function discard
 * current event information.
 *
 * NOTE: between ced_new_event() and ced_draw_event()
 *       must be some element creation calls.
 */
void ced_draw_event(void);

/*
 * This function really attempt to display event in CED.
 * Unlike ced_draw_event() does not reset the event.
 *
 * NOTE: between ced_new_event() and ced_draw_event()
 *       must be some element creation calls.
 */
void ced_send_event(void);

int ced_selected_id(void);

//hauke
int ced_selected_id_noblock(void);


/*********************************************
 *
 * The following is elements API.
 *
 *********************************************/

void ced_register_elements(void);

typedef enum {
  CED_TYPE_SHIFT=0x0,
  CED_LAYER_SHIFT=0x8
} CED_TYPE_BITS;

/*
typedef enum {
  CED_TYPE_SHIFT=0x0,
  CED_LAYER_SHIFT=0x0
} CED_TYPE_BITS;
*/


typedef struct {
  float x;
  float y;
  float z;
} CED_Point;

//class CED_Point{
//  public:
//  CED_Point(float _x, float _y, float _z){
//       x=_x;
//       y=_y;
//       z=_z;
//   }
//
//  CED_Point(void){
//  }
//  float x;
//  float y;
//  float z;
//}; 


/*
 * Hit element
 */

typedef enum {
    CED_HIT_POINT=0,
    CED_HIT_CROSS,
    CED_HIT_STAR,
    CED_HIT_BOX,
    CED_HIT_VXD
} CED_HIT_TYPE;

typedef struct {
  CED_Point p;
  unsigned type;  // point, star, etc
  unsigned layer; //layer
  unsigned color; // in ARGB form (so, 0xff0000 is RED)
  unsigned size;  // size of point/size of cross
  unsigned lcioID; // unique id of LICO object
} CED_Hit;

void ced_hit(float x,float y,float z,unsigned type,unsigned size,unsigned color);

//to give a bit of downward compatibility
void ced_hit_ID_old(float x,float y,float z,unsigned type, unsigned size,unsigned color, unsigned lcioID);


void ced_hit_ID(float x,float y,float z,unsigned type,unsigned layer, unsigned size,unsigned color, unsigned lcioID);

/*
 * Line element
 */

typedef struct {
  CED_Point p0;
  CED_Point p1;
  unsigned type;  // not yet defined...
  unsigned width; // not yet defined...
  unsigned color; // in ARGB form (so, 0xff0000 is RED)
  unsigned lcioID; // unique id of LICO object

} CED_Line;

void ced_line(float x0,float y0,float z0,
	      float x1,float y1,float z1,
	      unsigned type,unsigned width,unsigned color);
void ced_line_ID(float x0,float y0,float z0,
	      float x1,float y1,float z1,
	      unsigned type,unsigned width,unsigned color, unsigned lcioID);

/*
 * GeoCylinder
 */
typedef struct {
  float d;       // radius
  unsigned  sides;   // poligon order
  float rotate;  // angle degree
  float z;       // 1/2 length
  float shift;   // in z
  unsigned color;
} CED_GeoCylinder;

/*
 * GeoTube
 */
typedef struct {
  float r_o;            // outer radius
  float r_i;            // inner radius
  unsigned edges_o;     // edges outer
  unsigned edges_i;     // edges inner
  float rotate_o;       // angle degree, rotate outer cylinder
  float rotate_i;       //rotate inner cylinder
  float z;              // 1/2 length
  float shift;          // shift in z
  unsigned color;       // color
  unsigned type;        //describes the layer where this element lies
  bool classic_inner;   //draw the outer detector lines in classic view?
  bool classic_outer;   //draw the inner detector lines in classic view?
}  CED_GeoTube;


/** Same as CED_GeoTube but here as C++ struct with contstructor. This is allows
 *  to dynamically allocate the detector structure (in an std::vector using the constructor) 
 *  without knowing the exact number of  detector elements a priori (such as the #layers in the FTD).
 *  This wasn't poassible with the static allocation using the C type struct.
 *  ( Used in MArlinCED::drawGearDetector ).
 */
struct CEDGeoTube{
  CEDGeoTube(double  r_out,
	     double  r_in,
	     int edges_out,
	     int edges_in,
	     double  rotate_out,
	     double  rotate_in,
	     double  zlength,
	     double  zshift,
	     int col,
	     int layer,
	     bool classic_i,
	     bool classic_o ) :
    r_o(r_out),
    r_i(r_in),
    edges_o(edges_out),
    edges_i (edges_in),
    rotate_o (rotate_out),
    rotate_i (rotate_in),
    z( zlength),
    shift (zshift),
    color (col),
    type(layer),
    classic_inner(classic_i),
    classic_outer(classic_o) {} 
  float r_o;            // outer radius
  float r_i;            // inner radius
  unsigned edges_o;     // edges outer
  unsigned edges_i;     // edges inner
  float rotate_o;       // angle degree, rotate outer cylinder
  float rotate_i;       //rotate inner cylinder
  float z;              // 1/2 length
  float shift;          // shift in z
  unsigned color;       // color
  unsigned type;        //describes the layer where this element lies
  bool classic_inner;   //draw the outer detector lines in classic view?
  bool classic_outer;   //draw the inner detector lines in classic view?
} ;

void ced_geotubes(unsigned n,CED_GeoTube *all);


void ced_geocylinder(float d,unsigned sides,float rotate,float z,float shift,
		     unsigned color);

void ced_geocylinders(unsigned n,CED_GeoCylinder *all);
		   
/*
 * GeoCylinder rotatable
 * @author: S.Daraszewicz (UoE)
 * @date: 01.09.09
 */
typedef struct {
  float d;       	// radius
  unsigned sides; 	// poligon order
  float center[3];  // cylinder centre z,y,z
  float rotate[3];  // rotation angles wrt x,y,z axis
  float z;       	// length
  unsigned color;	// colour
  unsigned layer; 	// layer the Cylinder to be displayed onto
} CED_GeoCylinderR;
		   
		     
void ced_geocylinder_r(float d, double z, double * center, double * rotate, unsigned sides, 
		     unsigned int color, int layer);
		     

  /** GeoBox structure
   */
  typedef struct {
    /** The three box sizes in mm */
    double sizes[3];
    /** position of the center of the box*/
    double center[3];
    /** box color */
    unsigned int color;
  } CED_GeoBox;

  /** Send/Draw a box at position center (x,y,z in mm) with lengths along the 
   * axes specified in sizes.
   * 
   * @author A.Bulgheroni, INFN
   */
  void ced_geobox(double * sizes, double * center, unsigned int color );
  void ced_geobox_ID(double *size, double *position, unsigned int layer, unsigned int color, unsigned int lcio_id);
  void ced_geobox_r_ID(double *size, double *position, double *rotate, unsigned int layer, unsigned int color, unsigned int lcio_id);
  void rotate3d(double *vektor, double *rotate);



   /* 
   * @author A.Bulgheroni, INFN
   */
  void ced_geoboxes( unsigned int nBox, CED_GeoBox * allBoxes);

  typedef struct {
    /** The three box sizes in mm */
    double sizes[3];
    /** position of the center of the box*/
    double center[3];
    /** box color */
    unsigned int color;
    /** rotation angle in degrees */
    double rotate[3];
    /** layer for toggling display */
    unsigned int layer;
  } CED_GeoBoxR;

void ced_geobox_r(double * sizes, double * center, double * rotate, unsigned int color, unsigned int layer);
void ced_geobox_r_solid(double * sizes, double * center, double * rotate, unsigned int color, unsigned int layer);

//hauke
  typedef struct{
    char text[1000];
    int id;
  } CED_PICKING_TEXT; 

void ced_picking_text(const char *, int number);

//--------

  typedef struct{
    char text[400];
    int id;
  } CED_TEXT; 


void ced_describe_layer(const char *, int); //, int, int);


  typedef struct{
    char str[400];
    unsigned int id;
  } LAYER_TEXT; 

void ced_layer_text(char *, int);
//end hauke


/*
 * Energy spectrum colour map legend.
 * @author: S.Daraszewicz (UoE)
 * @date: 01.09.09
 */
  typedef struct {  
  	/** min energy on the legend */	
  	float ene_max;
  	/** max energy on the legend */
  	float ene_min;
  	/** number of ticks on the legend */
  	unsigned int ticks;
  	/** spectrum colour steps */
  	unsigned int color_steps; 
  	/** spectrum colour matrix */
  	unsigned int rgb_matrix[512][3]; //FIX ME: 512 size not changed with color_steps
  	/** LOG or LIN */
  	char scale;
  } CED_Legend;

void ced_legend(float ene_min, float ene_max, unsigned int color_steps, unsigned int ** rgb_matrix, unsigned int ticks, char scale);

  typedef struct {  
  	/** position of the centre of the base */	
  	double center[3];
  	/** rotation matrix */
  	double rotate[3];
    /** layer for toggling display */
    unsigned int layer;
    /** base radius */
    float base;
    /** height */
    float height;
    /** RGBA color */
    float RGBAcolor[4];
    unsigned lcioid; //hauke
  } CED_ConeR;


void ced_cone_r(float base, float height, double *center, double *rotate, unsigned int layer, float *RGBAcolor);
void ced_cone_r_ID(float base, float height, double *center, double *rotate, unsigned int layer, float *RGBAcolor, int lcioid); //hauke


  typedef struct {  
  	/** position of the centre of the base */	
  	double center[3];
  	/** rotation matrix */
  	double rotate[3];
    /** layer for toggling display */
    unsigned int layer;
    /** xyz size */
	double size[3];
    /** RGBA color */
   	int color;
    unsigned lcioid; //hauke
  } CED_EllipsoidR;


void ced_ellipsoid_r(double *size, double *center, double *rotate, unsigned int layer, int color);
void ced_ellipsoid_r_ID(double *size, double *center, double *rotate, unsigned int layer, int color, int lcioid); //hauke


  typedef struct {  
  	/** position of the centre of the base */	
  	double center[3];
  	/** rotation matrix */
  	double rotate[3];
    /** layer for toggling display */
    unsigned int layer;
    /** base radius */
	float radius;
	/** half height */
	float height;
    /** RGBA color */
    int color;
    unsigned lcioid; //hauke
  } CED_CluEllipseR;


void ced_cluellipse_r(float radius, float height, float *center, double *rotate, unsigned int layer, int color);
void ced_cluellipse_r_ID(float radius, float height, float *center, double *rotate, unsigned int layer, int color, int lcioid); //hauke


#ifdef __cplusplus
 }
#endif
	

#endif /* __CED_CLI_H */
