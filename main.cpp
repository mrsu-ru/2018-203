#include <cstdlib>
#include <cstring>
#include "lab.h"
#include <iostream>
#include "ivanovii.h"
#include "zhalninrv.h"
#include "scherbakovdv.h"
#include "polyakovda.h"
#include "bagapovar.h"
#include "tarasovams.h"
#include "serguninaes.h"
#include "syusinaev.h"
#include "fedyanovaam.h"
#include "borisovrs.h"
#include "nefedovms.h"
#include "kuznetsovais.h"
#include "grishaevaoov.h"
#include "bulychevaoa.h"
#include "salinaa.h"
#include "mescheryakovam.h"

void print_usage(char* name);


int main(int argc, char** argv)
{
  if (argc < 3) {
    print_usage(argv[0]);
    return 0;
  }

  lab *l = NULL;
  if (strcmp(argv[1], "salinaa") == 0) {
    l = new salinaa();
  }
  else if (strcmp(argv[1], "zhalninrv") == 0) {
    l = new zhalninrv();
  }
  else if (strcmp(argv[1],"nefedovms")==0) {
    l = new nefedovms();
  }
  else if (strcmp(argv[1],"bagapovar")==0) {
    l = new bagapovar();
  }
  else if (strcmp(argv[1],"borisovrs")==0) {
    l = new borisovrs();
  }
  else if (strcmp(argv[1],"fedyanovaam")==0) {
    l = new fedyanovaam();
  }
  else if (strcmp(argv[1],"kuznetsovais")==0) {
    l = new kuznetsovais();
  }
  else if (strcmp(argv[1],"serguninaes")==0) {
    l = new serguninaes();
  }
  else if (strcmp(argv[1],"syusinaev")==0) {
    l = new syusinaev();
  }
  else if (strcmp(argv[1],"tarasovams")==0) {
    l = new tarasovams();
  }
  else if (strcmp(argv[1],"scherbakovdv")==0) {
    l = new scherbakovdv();
  }
  else  if (strcmp(argv[1], "polyakovda") == 0) {
    l = new polyakovda();
  }
  else  if (strcmp(argv[1], "grishaevaov") == 0) {
    l = new grishaevaov();
  }
  else  if (strcmp(argv[1], "bulychevaoa") == 0) {
    l = new bulychevaoa();
  }
  else  if (strcmp(argv[1], "mescheryakovam") == 0) {
    l = new mescheryakovam();
  }
  else  {
    print_usage(argv[0]);
    return 0;
  }

  l->read_file();
  l->run(atoi(argv[2]));
  l->write_result();
  l->check_result();

  //delete l; // TODO:
  return 0;
}


void print_usage(char* name)
{
  std::cout << "Usage:\n\n  " << name << " <fio> <lab_number>\n";
}
