.intel_syntax noprefix

.set pbytes,32
.set plimbs,4
.global p
.global _p
p: _p:
    .quad 0x50ca4291ffffffff, 0xd1b004f94a5952c9, 0xb25c76a437728f3b, 0xa3091565b678a990

.inv_min_p_mod_r: /* -p^-1 mod 2^64 */
    .quad 0x50ca429200000001

.global fp_0
.global _fp_0
fp_0: _fp_0:
    .zero pbytes

.global fp_1
.global _fp_1
fp_1: _fp_1: /* 2^256 mod p */
    .quad 0xaf35bd6e00000001, 0x2e4ffb06b5a6ad36, 0x4da3895bc88d70c4, 0x5cf6ea9a4987566f

.r_squared_mod_p: /* (2^256)^2 mod p */
    .quad 0xce1eebd1f8f7c152, 0x3d0e272c4d770926, 0x868e6f58f9294cfb, 0x578b875020f84a87

.p_minus_2:
    .quad 0x50ca4291fffffffd, 0xd1b004f94a5952c9, 0xb25c76a437728f3b, 0xa3091565b678a990

.p_minus_1_halves:
    .quad 0xa8652148ffffffff, 0xe8d8027ca52ca964, 0x592e3b521bb9479d, 0x51848ab2db3c54c8

.p_plus_1_quarter:
	.quad 0x543290a480000000, 0xf46c013e529654b2, 0x2c971da90ddca3ce, 0x28c245596d9e2a64

.data
.global fp_mul_count
.global _fp_mul_count
fp_mul_count: _fp_mul_count:
    .quad 0

.text
.p2align 4,,15

.global fp_copy
.global _fp_copy
fp_copy: _fp_copy:
    cld
    mov rcx, plimbs
    rep movsq
    ret

.global fp_set
.global _fp_set
fp_set: _fp_set:
    push rdi
    call uintbig_set
    pop rdi
    mov rsi, rdi
    jmp fp_enc   
    
.global fp_add2
.global _fp_add2
fp_add2: _fp_add2:
  mov    rdx, rdi
.global fp_add3
.global _fp_add3
fp_add3: _fp_add3:
  push   r12  
  push   r15
  xor    rax, rax
  xor    r12, r12
  mov    r8, [rsi]
  mov    r9, [rsi+8]
  mov    r10, [rsi+16]
  mov    r11, [rsi+24]
  add    r8, [rdx] 
  adc    r9, [rdx+8] 
  adc    r10, [rdx+16] 
  adc    r11, [rdx+24]  
  adc    r12, 0
  mov    r15, [rip+p]
  sub    r8, r15
  mov    rcx, [rip+p+8]
  sbb    r9, rcx
  mov    rsi, [rip+p+16]
  sbb    r10, rsi
  mov    rdx, [rip+p+24]
  sbb    r11, rdx
  sbb    r12, 0
  sbb    rax, 0
  
  and    r15, rax
  and    rcx, rax
  and    rsi, rax
  and    rdx, rax
  
  add    r8, r15  
  adc    r9, rcx  
  adc    r10, rsi  
  adc    r11, rdx 
  mov    [rdi], r8
  mov    [rdi+8], r9 
  mov    [rdi+16], r10 
  mov    [rdi+24], r11
  pop    r15
  pop    r12
  ret
  
.global fp_sub2
.global _fp_sub2
fp_sub2: _fp_sub2:
  mov    rdx, rdi
  xchg   rsi, rdx
.global fp_sub3
.global _fp_sub3
fp_sub3: _fp_sub3:
  push   r12  
  xor    rax, rax
  mov    r8, [rsi]
  mov    r9, [rsi+8]
  mov    r10, [rsi+16]
  mov    r11, [rsi+24]
  sub    r8, [rdx] 
  sbb    r9, [rdx+8] 
  sbb    r10, [rdx+16] 
  sbb    r11, [rdx+24]
  sbb    rax, 0
  
  mov    r12, [rip+p]
  mov    rcx, [rip+p+8]
  mov    rsi, [rip+p+16]
  mov    rdx, [rip+p+24]
  and    r12, rax
  and    rcx, rax
  and    rsi, rax
  and    rdx, rax  
  add    r8, r12  
  adc    r9, rcx 
  adc    r10, rsi  
  adc    r11, rdx 
  mov    [rdi], r8
  mov    [rdi+8], r9 
  mov    [rdi+16], r10 
  mov    [rdi+24], r11 
  pop    r12
  ret

///////////////////////////////////////////////////////////////// MACROS
// z = a x bi + z
// Inputs: base memory pointer M1 (a),
//         bi pre-stored in rdx,
//         accumulator z in [Z0:Z4]
// Output: [Z0:Z4]
// Temps:  regs T0:T1
/////////////////////////////////////////////////////////////////
.macro MULADD64x256 M1, Z0, Z1, Z2, Z3, Z4, Z5, T0, T1, C, CC
    mulx   \T0, \T1, \M1     // A0*B0
    xor    \C, \C
    xor    \CC, \CC
    adox   \Z0, \T1
    adox   \Z1, \T0  
    mulx   \T0, \T1, 8\M1    // A0*B1
    adcx   \Z1, \T1
    adox   \Z2, \T0    
    mulx   \T0, \T1, 16\M1   // A0*B2
    adcx   \Z2, \T1
    adox   \Z3, \T0
    mulx   \T0, \T1, 24\M1   // A0*B3          
    adcx   \Z3, \T1
    adox   \Z4, \T0
    adox   \Z5, \CC
    adcx   \Z4, \CC 
    adc    \Z5, 0  
.endm
  
//***********************************************************************
//  Multiplication in GF(p^2), non-complex part
//  Operation: c [rdi] = a0 x b0 - a1 x b1
//  Inputs: a = [a1, a0] stored in [rsi] 
//          b = [b1, b0] stored in [rdx] 
//  Output: c stored in [rdi]
//***********************************************************************
.global fp2_mul_c0
fp2_mul_c0:    
    push   r12 
    push   r13 
    push   r14 
    push   r15
    mov    rcx, rdx
	
	// [rdi0:3] <- p - b1
	mov    r8, [rip+p]  
	mov    r9, [rip+p+8]   
	mov    r10, [rip+p+16]
	mov    r11, [rip+p+24] 
	mov    rax, [rcx+32]
	mov    rdx, [rcx+40]        
	sub    r8, rax
	sbb    r9, rdx
	mov    rax, [rcx+48]
	mov    rdx, [rcx+56]
	sbb    r10, rax
	sbb    r11, rdx
	mov    [rdi], r8
	mov    [rdi+8], r9
	mov    [rdi+16], r10
	mov    [rdi+24], r11
    
    // [r8:r12] <- z = a0 x b00 - a1 x b10
    mov    rdx, [rcx]
    mulx   r9, r8, [rsi]         
    xor    r13, r13
    mulx   r10, r11, [rsi+8]
    adox   r9, r11        
    mulx   r11, r12, [rsi+16] 
    adox   r10, r12        
    mulx   r12, r14, [rsi+24]
    adox   r11, r14  
    adox   r12, r13
           
    mov    rdx, [rdi]    
    MULADD64x256 [rsi+32], r8, r9, r10, r11, r12, r13, r14, r15, rax, rax 
    // [r9:r13] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r8             // rdx <- z0
    MULADD64x256 [rip+p], r8, r9, r10, r11, r12, r13, r14, r15, rax, rax
    
    // [r9:r13, r8] <- z = a0 x b01 - a1 x b11 + z 
    mov    rdx, [rcx+8]
    MULADD64x256 [rsi], r9, r10, r11, r12, r13, r8, r14, r15, r8, rax
    mov    rdx, [rdi+8]    
    MULADD64x256 [rsi+32], r9, r10, r11, r12, r13, r8, r14, r15, rax, rax
    // [r10:r13, r8] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r9             // rdx <- z0
    MULADD64x256 [rip+p], r9, r10, r11, r12, r13, r8, r14, r15, rax, rax
    
    // [r10:r13, r8:r9] <- z = a0 x b02 - a1 x b12 + z 
    mov    rdx, [rcx+16]
    MULADD64x256 [rsi], r10, r11, r12, r13, r8, r9, r14, r15, r9, rax
    mov    rdx, [rdi+16]    
    MULADD64x256 [rsi+32], r10, r11, r12, r13, r8, r9, r14, r15, rax, rax
    // [r11:r13, r8:r9] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r10            // rdx <- z0
    MULADD64x256 [rip+p], r10, r11, r12, r13, r8, r9, r14, r15, rax, rax
    
    // [r11:r13, r8:r10] <- z = a0 x b03 - a1 x b13 + z 
    mov    rdx, [rcx+24]
    MULADD64x256 [rsi], r11, r12, r13, r8, r9, r10, r14, r15, r10, rax
    mov    rdx, [rdi+24]    
    MULADD64x256 [rsi+32], r11, r12, r13, r8, r9, r10, r14, r15, rax, rax
    // [r12:r13, r8:r10] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r11            // rdx <- z0
    MULADD64x256 [rip+p], r11, r12, r13, r8, r9, r10, r14, r15, rax, rax

	// Final correction                        
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
    adc    r10, 0
    
    xor    rax, rax
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
 
    mov    [rdi], r12          
    mov    [rdi+8], r13         
    mov    [rdi+16], r8         
    mov    [rdi+24], r9  
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret
  
//***********************************************************************
//  Multiplication in GF(p^2), complex part
//  Operation: c [rdi] = a0 x b1 + a1 x b0
//  Inputs: a = [a1, a0] stored in [rsi] 
//          b = [b1, b0] stored in [rdx] 
//  Output: c stored in [rdi]
//***********************************************************************
.global fp2_mul_c1
fp2_mul_c1:    
    push   r12 
    push   r13 
    push   r14 
    push   r15  
    mov    rcx, rdx
    
    // [r8:r12] <- z = a0 x b10 + a1 x b00
    mov    rdx, [rcx+32]
    mulx   r9, r8, [rsi]         
    xor    r13, r13
    mulx   r10, r11, [rsi+8]
    adox   r9, r11        
    mulx   r11, r12, [rsi+16] 
    adox   r10, r12        
    mulx   r12, r14, [rsi+24]
    adox   r11, r14  
    adox   r12, r13
           
    mov    rdx, [rcx]    
    MULADD64x256 [rsi+32], r8, r9, r10, r11, r12, r13, r14, r15, rax, rax
    // [r9:r13] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r8             // rdx <- z0
    MULADD64x256 [rip+p], r8, r9, r10, r11, r12, r13, r14, r15, rax, rax
    
    // [r9:r13, r8] <- z = a0 x b01 - a1 x b11 + z 
    mov    rdx, [rcx+40]
    MULADD64x256 [rsi], r9, r10, r11, r12, r13, r8, r14, r15, r8, rax
    mov    rdx, [rcx+8]    
    MULADD64x256 [rsi+32], r9, r10, r11, r12, r13, r8, r14, r15, rax, rax
    // [r10:r13, r8] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r9             // rdx <- z0
    MULADD64x256 [rip+p], r9, r10, r11, r12, r13, r8, r14, r15, rax, rax
    
    // [r10:r13, r8:r9] <- z = a0 x b02 - a1 x b12 + z 
    mov    rdx, [rcx+48]
    MULADD64x256 [rsi], r10, r11, r12, r13, r8, r9, r14, r15, r9, rax
    mov    rdx, [rcx+16]    
    MULADD64x256 [rsi+32], r10, r11, r12, r13, r8, r9, r14, r15, rax, rax
    // [r11:r13, r8:r9] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r10            // rdx <- z0
    MULADD64x256 [rip+p], r10, r11, r12, r13, r8, r9, r14, r15, rax, rax
    
    // [r11:r13, r8:r10] <- z = a0 x b03 - a1 x b13 + z 
    mov    rdx, [rcx+56]
    MULADD64x256 [rsi], r11, r12, r13, r8, r9, r10, r14, r15, r10, rax
    mov    rdx, [rcx+24]    
    MULADD64x256 [rsi+32], r11, r12, r13, r8, r9, r10, r14, r15, rax, rax
    // [r12:r13, r8:r10] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   rax, rdx, r11            // rdx <- z0
    MULADD64x256 [rip+p], r11, r12, r13, r8, r9, r10, r14, r15, rax, rax

	// Final correction                        
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
    adc    r10, 0
    
    xor    rax, rax
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
   
    mov    [rdi], r12          
    mov    [rdi+8], r13         
    mov    [rdi+16], r8         
    mov    [rdi+24], r9  
    pop    r15
    pop    r14
    pop    r13
    pop    r12
    ret

///////////////////////////////////////////////////////////////// MACRO
// z = a x b (mod p)
// Inputs: base memory pointers M0 (a), M1 (b)
//         bi pre-stored in rdx,
//         accumulator z in [Z0:Z4], pre-stores a0 x b
// Output: [Z0:Z4]
// Temps:  regs T0:T1
/////////////////////////////////////////////////////////////////
.macro FPMUL256x256 M0, M1, Z0, Z1, Z2, Z3, Z4, Z5, T0, T1, C           
    // [Z1:Z5] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   \C, rdx, \Z0             // rdx <- z0
    MULADD64x256 [rip+p], \Z0, \Z1, \Z2, \Z3, \Z4, \Z5, \T0, \T1, \C, \C
    
    // [Z1:Z5, Z0] <- z = a01 x a1 + z 
    mov    rdx, 8\M0
    MULADD64x256 \M1, \Z1, \Z2, \Z3, \Z4, \Z5, \Z0, \T0, \T1, \Z0, \C
    // [Z2:Z5, Z0] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   \C, rdx, \Z1             // rdx <- z1
    MULADD64x256 [rip+p], \Z1, \Z2, \Z3, \Z4, \Z5, \Z0, \T0, \T1, \C, \C
    
    // [Z2:Z5, Z0:Z1] <- z = a02 x a1 + z  
    mov    rdx, 16\M0
    MULADD64x256 \M1, \Z2, \Z3, \Z4, \Z5, \Z0, \Z1, \T0, \T1, \Z1, \C
    // [Z3:Z5, Z0:Z1] <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   \C, rdx, \Z2             // rdx <- z2
    MULADD64x256 [rip+p], \Z2, \Z3, \Z4, \Z5, \Z0, \Z1, \T0, \T1, \C, \C
    
    // [Z3:Z5, Z0:Z2] <- z = a03 x a1 + z
    mov    rdx, 24\M0
    MULADD64x256 \M1, \Z3, \Z4, \Z5, \Z0, \Z1, \Z2, \T0, \T1, \Z2, \C
    // [Z4:Z5, Z0:Z2]  <- z = ((z0 x inv_min_p_mod_r mod 2^64) x p + z)/2^64
	mov    rdx, [rip+.inv_min_p_mod_r]
	mulx   \C, rdx, \Z3             // rdx <- z3
    MULADD64x256 [rip+p], \Z3, \Z4, \Z5, \Z0, \Z1, \Z2, \T0, \T1, \C, \C
.endm

//***********************************************************************
//  Squaring in GF(p^2), non-complex part
//  Operation: c [rdi] = (a0+a1) x (a0-a1)
//  Inputs: a = [a1, a0] stored in [rsi] 
//  Output: c stored in [rdi]
//***********************************************************************
.global fp2_sq_c0
fp2_sq_c0:   
    push   r12 
    push   r13
    push   r14

	// a0 + a1
    xor    r12, r12
	mov    rdx, [rsi]
	mov    r9, [rsi+8]
	mov    r10, [rsi+16]
	mov    r11, [rsi+24]
	add    rdx, [rsi+32]
	adc    r9, [rsi+40]
	adc    r10, [rsi+48]
	adc    r11, [rsi+56]
    adc    r12, 0

	// Correction
    xor    rax, rax
	mov    r8, [rip+p]
	mov    r13, [rip+p+8]
	mov    r14, [rip+p+16]
	mov    rcx, [rip+p+24]
	sub    rdx, r8
	sbb    r9, r13
	sbb    r10, r14
	sbb    r11, rcx
	sbb    r12, 0
	sbb    rax, 0
	and    r8, rax
	and    r13, rax
	and    r14, rax
	and    rcx, rax
	add    rdx, r8
	adc    r9, r13
	adc    r10, r14
	adc    r11, rcx
	mov    [rdi], rdx
	mov    [rdi+8], r9
	mov    [rdi+16], r10
	mov    [rdi+24], r11 
	
	// a0 - a1 + p
    xor    r14, r14
	mov    r8, [rsi]
	mov    r10, [rsi+8]
	mov    r12, [rsi+16]
	mov    r13, [rsi+24]
	sub    r8, [rsi+32]
	sbb    r10, [rsi+40]
	sbb    r12, [rsi+48] 
	sbb    r13, [rsi+56]
    sbb    r14, 0
	add    r8, [rip+p]                    
	adc    r10, [rip+p+8]
	adc    r12, [rip+p+16]
	adc    r13, [rip+p+24]
    adc    r14, 0

	// Correction
    xor    rax, rax
	mov    r9, [rip+p]
	mov    r11, [rip+p+8]
	mov    rcx, [rip+p+16]
	mov    rdx, [rip+p+24]
	sub    r8, r9
	sbb    r10, r11
	sbb    r12, rcx
	sbb    r13, rdx
	sbb    r14, 0
	sbb    rax, 0
	and    r9, rax
	and    r11, rax
	and    rcx, rax
	and    rdx, rax
	add    r8, r9
	adc    r10, r11
	adc    r12, rcx
	adc    r13, rdx
	mov    [rdi+32], r8               
	mov    [rdi+40], r10 
	mov    [rdi+48], r12 
	mov    [rdi+56], r13 
    
    // [r8:r12] <- z = a00 x a1
    mov    rdx, [rdi]
    mulx   r9, r8, r8        
    xor    rax, rax
    mulx   r10, r11, r10  
    adox   r9, r11        
    mulx   r11, r12, r12  
    adox   r10, r12        
    mulx   r12, r14, r13 
    mov    r13, 0
    adox   r11, r14
    adox   r12, rax 

    FPMUL256x256 [rdi], [rdi+32], r8, r9, r10, r11, r12, r13, r14, rcx, rax

	// Final correction                        
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
   
    mov    [rdi], r12          
    mov    [rdi+8], r13         
    mov    [rdi+16], r8         
    mov    [rdi+24], r9  
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Squaring in GF(p^2), complex part
//  Operation: c [rdi] = 2a0 x a1
//  Inputs: a = [a1, a0] stored in [reg_p1] 
//  Output: c stored in [rdi]
//***********************************************************************
.global fp2_sq_c1
fp2_sq_c1:
    push   r12 
    push   r13
    push   r14

	// 2*a0
    xor    r12, r12
	mov    rdx, [rsi]
	mov    r9, [rsi+8]
	mov    r10, [rsi+16]
	mov    r11, [rsi+24]
	add    rdx, rdx
	adc    r9, r9
	adc    r10, r10
	adc    r11, r11
    adc    r12, 0

	// Correction
    xor    rax, rax
	mov    r8, [rip+p]
	mov    r13, [rip+p+8]
	mov    r14, [rip+p+16]
	mov    rcx, [rip+p+24]
	sub    rdx, r8
	sbb    r9, r13
	sbb    r10, r14
	sbb    r11, rcx
	sbb    r12, 0
	sbb    rax, 0
	and    r8, rax
	and    r13, rax
	and    r14, rax
	and    rcx, rax
	add    rdx, r8
	adc    r9, r13
	adc    r10, r14
	adc    r11, rcx
	sub    rsp, 32
	mov    [rsp+8], r9
	mov    [rsp+16], r10 
	mov    [rsp+24], r11     
    
    // [r8:r12] <- z = a00 x a1
    mulx   r9, r8, [rsi+32]
    xor    r13, r13 
    mulx   r10, r11, [rsi+40]
    adox   r9, r11        
    mulx   r11, r12, [rsi+48]
    adox   r10, r12        
    mulx   r12, r14, [rsi+56]
    adox   r11, r14  
    adox   r12, r13  

	FPMUL256x256 [rsp], [rsi+32], r8, r9, r10, r11, r12, r13, r14, rcx, rax
	add    rsp, 32

	// Final correction                        
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
   
    mov    [rdi], r12          
    mov    [rdi+8], r13         
    mov    [rdi+16], r8         
    mov    [rdi+24], r9  
    pop    r14
    pop    r13
    pop    r12
    ret

//***********************************************************************
//  Field multiplication in GF(p)
//  Operation: c = a x b mod p
//  Inputs: a stored in [rsi], b stored in [rdx] 
//  Output: c stored in [rdi]
//***********************************************************************
.global fp_mul2
.global _fp_mul2
fp_mul2: _fp_mul2:
    mov    rdx, rdi
.global fp_mul3
.global _fp_mul3
fp_mul3: _fp_mul3: 
    push   r12
    push   r13 
    push   r14 
    push   r15 
    mov    rcx, rdx 
     
    // [r8:r12] <- z = a x b0
    mov    rdx, [rcx]
    mulx   r9, r8, [rsi]
    xor    r13, r13 
    mulx   r10, r11, [rsi+8]
    adox   r9, r11        
    mulx   r11, r12, [rsi+16]
    adox   r10, r12        
    mulx   r12, r14, [rsi+24] 
    adox   r11, r14
    adox   r12, r13 

	FPMUL256x256 [rcx], [rsi], r8, r9, r10, r11, r12, r13, r14, r15, rax

	// Final correction                        
	mov    rsi, [rip+p]
	mov    rcx, [rip+p+8]
	mov    rdx, [rip+p+16]
	mov    r11, [rip+p+24]
	sub    r12, rsi
	sbb    r13, rcx
	sbb    r8, rdx
	sbb    r9, r11
	sbb    r10, 0
	sbb    rax, 0
	and    rsi, rax
	and    rcx, rax
	and    rdx, rax
	and    r11, rax
	add    r12, rsi
	adc    r13, rcx
	adc    r8, rdx
	adc    r9, r11
    
    mov    [rdi], r12          
    mov    [rdi+8], r13         
    mov    [rdi+16], r8         
    mov    [rdi+24], r9 
    pop    r15  
    pop    r14
    pop    r13
    pop    r12
    ret

.global fp_cswap
.global _fp_cswap
fp_cswap: _fp_cswap:
    movzx rax, dl
    neg rax
    .set k, 0
    .rept plimbs
        mov rcx, [rdi + 8*k]
        mov rdx, [rsi + 8*k]

        mov r8, rcx
        xor r8, rdx
        and r8, rax

        xor rcx, r8
        xor rdx, r8

        mov [rdi + 8*k], rcx
        mov [rsi + 8*k], rdx

        .set k, k+1
    .endr
    ret

.reduce_once:
    push rbp
    mov rbp, rdi

    mov rdi, [rbp +  0]
    sub rdi, [rip + p +  0]
    mov rsi, [rbp +  8]
    sbb rsi, [rip + p +  8]
    mov rdx, [rbp + 16]
    sbb rdx, [rip + p + 16]
    mov rcx, [rbp + 24]
    sbb rcx, [rip + p + 24]
    sbb rax, 0 /* handle carry from caller */

    setnc al
    movzx rax, al
    neg rax

.macro cswap2, r, m
    xor \r, \m
    and \r, rax
    xor \m, \r
.endm

    cswap2 rdi, [rbp +  0]
    cswap2 rsi, [rbp +  8]
    cswap2 rdx, [rbp + 16]
    cswap2 rcx, [rbp + 24]

    pop rbp
    ret

/* Montgomery arithmetic */

.global fp_enc
.global _fp_enc
fp_enc: _fp_enc:
    lea rdx, [rip + .r_squared_mod_p]
    jmp fp_mul3

.global fp_dec
.global _fp_dec
fp_dec: _fp_dec:
    lea rdx, [rip + uintbig_1]
    jmp fp_mul3
    
.global fp_mul2_old
.global _fp_mul2_old
fp_mul2_old: _fp_mul2_old:
    mov rdx, rdi
.global fp_mul3_old
.global _fp_mul3_old
fp_mul3_old: _fp_mul3_old:
    push rbp
    push rbx
    push r12
    push rdi

    /* inc qword ptr fp_mul_count */

    mov rdi, rsi
    mov rsi, rdx

    xor r8,  r8
    xor r9,  r9
    xor r10, r10
    xor r11, r11
    xor r12, r12
    xor rbp, rbp

    /* flags are already cleared */

.macro MULSTEP, k, r0, r1, r2, r3, r4, r5

    mov rdx, [rsi +  0]
    mulx rcx, rdx, [rdi + 8*\k]
    add rdx, \r0
    mulx rcx, rdx, [rip + .inv_min_p_mod_r]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rip + p +  0]
    adox \r0, rax

    mulx rcx, rax, [rip + p +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rip + p + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rip + p + 24]
    adcx \r3, rbx
    adox \r3, rax

    mov rax, 0
    adcx \r4, rcx
    adox \r4, rax

    adcx \r5, rax
    adox \r5, rax

    mov rdx, [rdi + 8*\k]

    xor rax, rax /* clear flags */

    mulx rbx, rax, [rsi +  0]
    adox \r0, rax

    mulx rcx, rax, [rsi +  8]
    adcx \r1, rbx
    adox \r1, rax

    mulx rbx, rax, [rsi + 16]
    adcx \r2, rcx
    adox \r2, rax

    mulx rcx, rax, [rsi + 24]
    adcx \r3, rbx
    adox \r3, rax

    mov rax, 0
    adcx \r4, rcx
    adox \r4, rax

    adcx \r5, rax
    adox \r5, rax

.endm

    MULSTEP 0, r8,  r9,  r10, r11, r12, rbp
    MULSTEP 1, r9,  r10, r11, r12, rbp, r8
    MULSTEP 2, r10, r11, r12, rbp, r8,  r9
    MULSTEP 3, r11, r12, rbp, r8,  r9,  r10

    pop rdi

    mov [rdi +  0], r12
    mov [rdi +  8], rbp
    mov [rdi + 16], r8
    mov [rdi + 24], r9
    mov rax, r10

    pop r12
    pop rbx
    pop rbp
    jmp .reduce_once

.global fp_sq1
.global _fp_sq1
fp_sq1: _fp_sq1:
    mov rsi, rdi
.global fp_sq2
.global _fp_sq2
fp_sq2: _fp_sq2:
    /* TODO implement optimized Montgomery squaring */
    mov rdx, rsi
    jmp fp_mul3

/* (obviously) not constant time in the exponent! */
.fp_pow:
    push rbx
    mov rbx, rsi
    push r12
    push r13
    push rdi
    sub rsp, pbytes

    mov rsi, rdi
    mov rdi, rsp
    call fp_copy

    mov rdi, [rsp + pbytes]
    lea rsi, [rip + fp_1]
    call fp_copy

.macro POWSTEP, k
        mov r13, [rbx + 8*\k]
        xor r12, r12

        2:
        test r13, 1
        jz 1f

        mov rdi, [rsp + pbytes]
        mov rsi, rsp
        call fp_mul2

        1:
        mov rdi, rsp
        call fp_sq1

        shr r13

        inc r12
        test r12, 64
        jz 2b
.endm

    POWSTEP 0
    POWSTEP 1
    POWSTEP 2
    POWSTEP 3

    add rsp, pbytes+8
    pop r13
    pop r12
    pop rbx
    ret

/* TODO use a better addition chain? */
.global fp_inv
.global _fp_inv
fp_inv: _fp_inv:
    lea rsi, [rip + .p_minus_2]
    jmp .fp_pow

.global fp_sqrt
.global _fp_sqrt
fp_sqrt: _fp_sqrt:
    lea rsi, [rip + .p_plus_1_quarter]
    jmp .fp_pow

/* TODO use a better addition chain? */
.global fp_issquare
.global _fp_issquare
fp_issquare: _fp_issquare:
    push rdi
    lea rsi, [rip + .p_minus_1_halves]
    call .fp_pow
    pop rdi

    xor rax, rax
    .set k, 0
    .rept plimbs
        mov rsi, [rdi + 8*k]
        xor rsi, [rip + fp_1 + 8*k]
        or rax, rsi
        .set k, k+1
    .endr
    test rax, rax
    setz al
    movzx rax, al
    ret

/* not constant time (but this shouldn't leak anything of importance) */
.global fp_random
.global _fp_random
fp_random: _fp_random:

    push rdi
    mov rsi, pbytes
    call _randombytes
    pop rdi

    .set k, plimbs-1
    .rept plimbs
        mov rax, [rip + p + 8*k]
        cmp [rdi + 8*k], rax
        ja fp_random
        jb 0f
        .set k, k-1
    .endr
    jmp fp_random
    0:
    ret
