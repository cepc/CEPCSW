/* "C" event display.
 * Communications related part. 
 *
*ik
 * Alexey Zhelezov, DESY/ITEP, 2005 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>

#include <ced.h>

//hauke
//#include <stropts.h>
#include <poll.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <netdb.h>
#include <sys/socket.h> /* for AF_INET */
#include <iostream>


//http://www.rhyolite.com/pipermail/dcc/2004/001986.html
#ifndef POLLRDNORM //fg: should be defined in poll.h
# define POLLRDNORM     0x040           /* Normal data may be read.  */
# define POLLRDBAND     0x080           /* Priority data may be read.  */
# define POLLWRNORM     0x100           /* Writing now will not block.  */
# define POLLWRBAND     0x200           /* Priority data may be written.  */
#endif
//end hauke

static int ced_fd=-1; // CED connection socket

static unsigned short ced_port=7927; // port No of CED (assume localhost)
static char ced_host[30];

// Return 0 if can be connected, -1 otherwise.
/*static*/ int ced_connect(void){
  static time_t last_attempt=0;
  time_t ct;
  struct sockaddr_in addr;

  if(ced_fd>=0){
    return 0; // already connected;
  }
  time(&ct);
  if(ct-last_attempt<5){
    return -1; // don't try reconnect all the time
  }
  addr.sin_family=AF_INET;
  addr.sin_port=htons(ced_port);
  addr.sin_addr.s_addr=inet_addr(ced_host); 
  memset(&addr.sin_zero, 0, sizeof(addr.sin_zero)); //not nessesary because sin_zero is not used!

  ced_fd=socket(PF_INET,SOCK_STREAM,0);
  if(connect(ced_fd,(struct sockaddr *)&addr,sizeof(addr)) != 0){
    if(!last_attempt){
        perror("WARNING:CED: can't connect to CED");
    }
    time(&last_attempt);
    close(ced_fd);
    ced_fd=-1;
    return -1;
  }
  fprintf(stderr,"INFO:CED: connected to CED\n");
  return 0;
}


typedef struct {
  unsigned size;            // size of one item in bytes
  unsigned char *b;         // "body" - data are stored here
                            // (here is some trick :)
  unsigned long count;      // number of usefull items
  unsigned long alloced;    // number of allocated items
  ced_draw_cb draw;         // draw fucation, NOT used in CED client
} ced_element;

typedef struct {
  ced_element *e;
  unsigned      e_count;
} ced_event;

//static ced_event eve = {0,0};

static ced_event eve = {0,0};

// NOT used in CED client
static ced_event ceve = {0,0}; // current event on screen

// we reserve this size just before ced_element.b data
#define HDR_SIZE 8 

unsigned ced_register_element(unsigned item_size,ced_draw_cb draw_func){
  ced_element *pe;
  if(!(eve.e_count&0xf)){
    eve.e=(ced_element *) realloc(eve.e,(eve.e_count+0x10)*sizeof(ced_element));
  }

  pe=eve.e+eve.e_count;
  memset(pe,0,sizeof(*pe));
  pe->size=item_size;
  pe->draw=draw_func;
  return eve.e_count++;
}

static void ced_reset(void){
  unsigned i;
  
  for(i=0;i<eve.e_count;i++){
    eve.e[i].count=0;
   // if( eve.e[i].alloced > 0){
   //     eve.e[i].alloced=0; //hauke: 15.12.11
   //     //free(eve.e[i].b-HDR_SIZE);
   // }
  }
}

static void ced_buf_alloc(ced_element *pe,unsigned count){
  if(!pe->b){

    //std::cout << "malloc  requestet: " << count*pe->size+HDR_SIZE << "bytes" << std::endl;
    pe->b=(unsigned char *) malloc(count*pe->size+HDR_SIZE);
    //printf("malloc: ask for NEW %lu bytes pointer: %p\n ", count*pe->size+HDR_SIZE, pe->b); //hauke
    if(pe->b==NULL){ //hauke
        printf("ERROR: malloc failed!\n");
        exit(1);
    }
  }else{
    //free(pe->b-HDR_SIZE);
    //pe->b=(unsigned char *) malloc(count*pe->size+HDR_SIZE);

    //std::cout << "realloc requestet: " << count*pe->size+HDR_SIZE << "bytes" << std::endl;
    pe->b=(unsigned char *) realloc(pe->b-HDR_SIZE,count*pe->size+HDR_SIZE);

    //printf("malloc: ask for %lu bytes, pointer: %p\n", count*pe->size+HDR_SIZE,pe->b);//hauke
    if(pe->b==NULL){ //hauke
        printf("ERROR: malloc failed!\n");
        exit(1);
    }
  }
  pe->b+=HDR_SIZE;
  pe->alloced=count;
}

void *ced_add(unsigned id){
  ced_element *pe;
  if(id >= eve.e_count){
    fprintf(stderr,"BUG:CED: attempt to access not registered element\n");
    return 0;
  }
  pe=eve.e+id;

  if(pe->count==pe->alloced){
    ced_buf_alloc(pe,pe->alloced+256);
  }

  return (pe->b+(pe->count++)*pe->size);
}

static void ced_event_copy(ced_event *trg){
  unsigned i;
  ced_element *pe;
  //std::cout << "trg->e_count: " << trg->e_count << std::endl;
  //std::cout << "eve.e_count: " << eve.e_count << std::endl;

  //eve.e_count = 0;
  if(trg->e_count<eve.e_count){
    //free(trg->e);
    trg->e=(ced_element*) realloc(trg->e,eve.e_count*sizeof(ced_element));

    //trg->e=(ced_element*) malloc(eve.e_count*sizeof(ced_element));
  }


  for(i=0;i<eve.e_count;i++){
        pe=trg->e+i;
        if(i<trg->e_count){
            //if(pe->alloced > 0){
            //  free(pe->b);
            //  pe->b=NULL;
            //  pe->alloced=0;
            //}


            //if(pe->alloced > 0){
            //    std::cout << "try to free" << std::endl;
            //    free(pe->b-HDR_SIZE);
            //    pe->alloced = 0;
            //    std::cout << "finished" << std::endl;
            //}

            if(pe->alloced<eve.e[i].alloced) {
	          ced_buf_alloc(pe,eve.e[i].alloced);
                //std::cout << "test1 " << std::endl;
            }
            pe->count=eve.e[i].count;

        }else{
            memcpy(pe,eve.e+i,sizeof(ced_element));
            if(pe->b){
	              pe->b=0;
                  //  std::cout << "test2 " << std::endl;
	              ced_buf_alloc(pe,pe->alloced);
            }
        }
        if(pe->count){
            memcpy(pe->b,eve.e[i].b,pe->count*pe->size);
        }
  }
  trg->e_count=eve.e_count;
}

void ced_do_draw_event(void){
  unsigned int i,j;
  ced_element *pe;
  unsigned char *pdata;
  for(i=0;i<ceve.e_count;i++){
    //printf("ceve.e_count: %i\n", ceve.e_count);
    //for(i=ceve.e_count-1; i >=0;i--){ //quick hack, change order so that the detector is drawn at last
    //printf("i = %i\n", i);

    pe=ceve.e+i;
    if(!pe->draw)
      continue;
    for(pdata=pe->b,j=0;j<pe->count;j++,pdata+=pe->size)
      (*(pe->draw))(pdata);
  }
}

typedef enum {
  DRAW_EVENT=10000
} MSG_TYPE;

int ced_process_input(void *data){
  struct _phdr{
    unsigned size;
    unsigned type;
    unsigned char b[4];
  } *hdr = (_phdr*) data;
  unsigned count;
  ced_element *pe;
  
  if(!data){ // new client is connected
    ced_reset();
    return 0;
  }

  if(hdr->type == DRAW_EVENT){
    ced_event_copy(&ceve);
    ced_reset();
    return 1;
  }
  if(hdr->type>=eve.e_count){
    fprintf(stderr,"WARNING:CED: undefined element type (%u), ignored\n",
	    hdr->type);
    return 0;
  }
  pe=eve.e+hdr->type;
  if((hdr->size-HDR_SIZE)%pe->size){
    fprintf(stderr,"BUG:CED: size alignment is wrong for element %u\n", hdr->type);
    return 0;
  }
  count=(hdr->size-HDR_SIZE)/pe->size;
  if(!count)
    return 0;
  if(count>=pe->alloced)
    ced_buf_alloc(pe,count+256);
  memcpy(pe->b,hdr->b,count*pe->size);
  pe->count=count;
  return 0;
}

void ced_send_event(void){
  struct _phdr{
    int size;
    unsigned type;
  } *hdr,draw_hdr;
  unsigned i,problem=0;
  int sent_sum;
  char *buf;
  int sent;
  ced_element *pe;

  if(ced_connect())
    return;
  for(i=0;i<eve.e_count && !problem;i++){
    //printf("i=%i\n",i);
    pe=eve.e+i;
    if(!pe->count)
      continue;
    
    //unsigned hauke;
    //printf("size of unsigned %i\n", sizeof(hauke));
    hdr=(struct _phdr *)(pe->b-HDR_SIZE); // !!! HERE is the trick :)
    hdr->type=i;
    //printf("pe->count %i, pe->size %i\n",pe->count, pe->size);
    hdr->size=HDR_SIZE+pe->count*pe->size;
    sent_sum=0;
    //if(hdr->size > 10000000){printf("U P S!  This data set is realy big! (%f kB)(%i counts)\n",(hdr->size)/1024.0,pe->count);}
    buf=(char *)hdr;
    //printf("hdr->size=%i\n",hdr->size);
    while(sent_sum<hdr->size){
        //printf("sent_sum = %i, hdr->size=%i\n",sent_sum,hdr->size);
	    sent=write(ced_fd,buf+sent_sum,hdr->size-sent_sum);
        
        //printf("byte: %u\n", buf[sent_sum]);
	    if(sent<0){
            printf("send < 0\n");
	        problem=1;
	        break;
	    }
	    sent_sum+=sent;
    }
  }
  if(!problem){
    draw_hdr.size=HDR_SIZE;
    draw_hdr.type=DRAW_EVENT;
    if(write(ced_fd,&draw_hdr,HDR_SIZE)!=HDR_SIZE)
      problem=1;
  }
  if(problem){
    perror("WARNING:CED: can't send event, till next time...");
    close(ced_fd);
    ced_fd=-1;
  }
}


//hauke
int ced_selected_id_noblock() {
  int id=-1 ;
  struct pollfd fds[1];
  fds[0].fd=ced_fd;
  fds[0].events = POLLRDNORM | POLLIN;
  if(poll(fds,1,0) > 0){
    if(recv(ced_fd, &id, sizeof(int) , 0 ) > 0){
        return id;
    }else{
        return -1;
    }
  }else{
   return -1;
  }
}

int ced_selected_id() {
  int id=-1 ;
  if(recv(ced_fd, &id, sizeof(int) , 0 ) > 0){
     return id;
  }else{
     return -1;
  }
}
#include <signal.h>
// API
void ced_client_init(const char *hostname,unsigned short port){
  struct hostent *host = gethostbyname(hostname);
  snprintf(ced_host, 30, "%u.%u.%u.%u\n",(unsigned char)host->h_addr[0] ,(unsigned char)host->h_addr[1] ,(unsigned char)host->h_addr[2] ,(unsigned char)host->h_addr[3]); 


  //printf("ip: %s\n",  ced_host);
  //ced_host=host->h_addr;
  ced_port=port;
  signal(SIGPIPE,SIG_IGN);
}

void ced_new_event(void){
  ced_reset();
}

void ced_draw_event(void){
  ced_send_event();
  ced_reset();
}
