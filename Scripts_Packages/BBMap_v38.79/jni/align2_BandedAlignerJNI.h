/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class align2_BandedAlignerJNI */

#ifndef _Included_align2_BandedAlignerJNI
#define _Included_align2_BandedAlignerJNI
#ifdef __cplusplus
extern "C" {
#endif
#undef align2_BandedAlignerJNI_big
#define align2_BandedAlignerJNI_big 99999999L
/*
 * Class:     align2_BandedAlignerJNI
 * Method:    alignForwardJNI
 * Signature: ([B[BIIIZI[B[I)I
 */
JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignForwardJNI
  (JNIEnv *, jobject, jbyteArray, jbyteArray, jint, jint, jint, jboolean, jint, jbyteArray, jintArray);

/*
 * Class:     align2_BandedAlignerJNI
 * Method:    alignForwardRCJNI
 * Signature: ([B[BIIIZI[B[B[I)I
 */
JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignForwardRCJNI
  (JNIEnv *, jobject, jbyteArray, jbyteArray, jint, jint, jint, jboolean, jint, jbyteArray, jbyteArray, jintArray);

/*
 * Class:     align2_BandedAlignerJNI
 * Method:    alignReverseJNI
 * Signature: ([B[BIIIZI[B[I)I
 */
JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignReverseJNI
  (JNIEnv *, jobject, jbyteArray, jbyteArray, jint, jint, jint, jboolean, jint, jbyteArray, jintArray);

/*
 * Class:     align2_BandedAlignerJNI
 * Method:    alignReverseRCJNI
 * Signature: ([B[BIIIZI[B[B[I)I
 */
JNIEXPORT jint JNICALL Java_align2_BandedAlignerJNI_alignReverseRCJNI
  (JNIEnv *, jobject, jbyteArray, jbyteArray, jint, jint, jint, jboolean, jint, jbyteArray, jbyteArray, jintArray);

#ifdef __cplusplus
}
#endif
#endif
