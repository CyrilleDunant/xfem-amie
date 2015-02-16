
#ifndef GL_UTILS_H
#define GL_UTILS_H

#include <GL/gl.h>
#include <cmath>

void glPerspective( GLdouble fovy, GLdouble aspect,
                              GLdouble zNear, GLdouble zFar ) ;

void glLookAt( GLdouble eyex, GLdouble eyey, GLdouble eyez,
                         GLdouble centerx, GLdouble centery, GLdouble
centerz,
                         GLdouble upx, GLdouble upy, GLdouble upz ) ;

#endif