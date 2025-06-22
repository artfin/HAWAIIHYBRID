#include <stdio.h>
#include "hawaii.h"

int main()
{
    String_Builder sb = {0};
   
    {
      printf("---------------------------------------------\n"); 
      sb_append(&sb, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", 257);
      printf("Checking that buffer sufficiently extends in 'sb_append' when the length of the appended string exceeds SB_INIT_CAPACITY:\n");
      printf("sb.items: '%s'\n", sb.items);
      printf("sb.count = %zu\n", sb.count);
      printf("sb.capacity = %zu\n", sb.capacity);
      printf("---------------------------------------------\n\n"); 

      sb_reset(&sb);
      sb_free(&sb);
    }

    {
      printf("---------------------------------------------\n"); 
      sb_append_cstring(&sb, "Hello, World");

      printf("Demonstrating that a C-style string is appended into String_Builder without terminating null-byte:\n");
      printf("sb.items: '%s'\n", sb.items);
      printf("sb.count = %zu\n", sb.count);
      printf("sb.capacity = %zu\n", sb.capacity);

      sb_append_cstring(&sb, "And also some string");
      printf("sb.items: '%s'\n", sb.items);

      sb_reset(&sb);
      printf("---------------------------------------------\n\n"); 
    }

    {
      printf("---------------------------------------------\n"); 
      printf("Checking that buffer sufficiently extends in 'sb_append_cstring' when the length of the appended string exceeds SB_INIT_CAPACITY:\n");
      sb_append_cstring(&sb, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");  
      printf("After appending very long string to String Builder:\n");
      printf("sb.items: '%s'\n", sb.items);
      printf("sb.count = %zu\n", sb.count);
      printf("sb.capacity = %zu\n", sb.capacity);
      printf("---------------------------------------------\n\n"); 
    }

    {
      printf("---------------------------------------------\n");
      sb_free(&sb);
      printf("After freeing String Builder:\n");
      printf("sb.items: '%s'\n", sb.items);
      printf("sb.count = %zu\n", sb.count);
      printf("sb.capacity = %zu\n", sb.capacity);
      printf("---------------------------------------------\n\n"); 
    }

    return 0;
}
