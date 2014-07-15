	.file	"prof_charge_assignment.c"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC5:
	.string	"Ik charge assignment took %lf s.\n"
	.align 8
.LC6:
	.string	"Ik force assignment took %lf s.\n"
	.align 8
.LC7:
	.string	"Ad charge assignment took %lf s.\n"
	.align 8
.LC8:
	.string	"Ad force assignment took %lf s.\n"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC9:
	.string	"Total time:"
.LC11:
	.string	"Ik: %lf\n"
.LC12:
	.string	"Ad: %lf\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB19:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movl	$10, %edx
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$112, %rsp
	.cfi_def_cfa_offset 144
	movq	8(%rsi), %rdi
	xorl	%esi, %esi
	call	strtol
	movsd	.LC1(%rip), %xmm1
	movl	%eax, %esi
	movsd	.LC0(%rip), %xmm0
	xorl	%edi, %edi
	call	generate_system
	movl	8(%rax), %edi
	movq	%rax, %rbx
	call	Init_forces
	movq	16(%rbp), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	movq	%rax, %r12
	call	strtol
	movq	24(%rbp), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	movl	%eax, 88(%rsp)
	call	strtol
	movl	%eax, %edx
	movl	%eax, 96(%rsp)
	leaq	64(%rsp), %rsi
	leal	-1(%rdx), %eax
	movq	%rbx, %rdi
	movl	%eax, 92(%rsp)
	movl	%edx, %eax
	imull	%edx, %eax
	imull	%edx, %eax
	movabsq	$4602678819172646912, %rdx
	movq	%rdx, 64(%rsp)
	movl	%eax, 100(%rsp)
	movabsq	$4613937818241073152, %rax
	movq	%rax, 72(%rsp)
	call	Init_ik
	movq	%rax, %rbp
	call	MPI_Wtime
	leaq	64(%rsp), %rsi
	xorl	%ecx, %ecx
	movq	%rbp, %rdx
	movq	%rbx, %rdi
	movsd	%xmm0, (%rsp)
	call	assign_charge
	call	MPI_Wtime
	movsd	(%rsp), %xmm1
	movl	$.LC5, %edi
	movl	$1, %eax
	subsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	xorpd	%xmm1, %xmm1
	movapd	%xmm0, %xmm2
	addsd	%xmm1, %xmm2
	movsd	%xmm1, (%rsp)
	movsd	%xmm2, 16(%rsp)
	call	printf
	call	MPI_Wtime
	leaq	64(%rsp), %rsi
	movapd	%xmm0, %xmm3
	movsd	.LC1(%rip), %xmm0
	xorl	%r8d, %r8d
	movq	%r12, %rcx
	movq	%rbp, %rdx
	movq	%rbx, %rdi
	movsd	%xmm3, 32(%rsp)
	call	assign_forces
	call	MPI_Wtime
	movsd	32(%rsp), %xmm3
	movl	$.LC6, %edi
	movsd	16(%rsp), %xmm2
	movl	$1, %eax
	subsd	%xmm0, %xmm3
	addsd	%xmm3, %xmm2
	movapd	%xmm3, %xmm0
	movsd	%xmm2, 16(%rsp)
	call	printf
	movq	%rbp, %rdi
	call	Free_data
	leaq	64(%rsp), %rsi
	movq	%rbx, %rdi
	call	Init_ad
	movq	%rax, %rbp
	call	MPI_Wtime
	leaq	64(%rsp), %rsi
	xorl	%ecx, %ecx
	movq	%rbp, %rdx
	movl	$1, %r8d
	movq	%rbx, %rdi
	movsd	%xmm0, 32(%rsp)
	call	assign_charge_and_derivatives
	call	MPI_Wtime
	movsd	32(%rsp), %xmm3
	movl	$.LC7, %edi
	movsd	(%rsp), %xmm1
	movl	$1, %eax
	subsd	%xmm0, %xmm3
	addsd	%xmm3, %xmm1
	movapd	%xmm3, %xmm0
	movsd	%xmm1, 56(%rsp)
	call	printf
	call	MPI_Wtime
	leaq	64(%rsp), %rsi
	movapd	%xmm0, %xmm1
	movsd	.LC1(%rip), %xmm0
	movq	%r12, %rcx
	movq	%rbp, %rdx
	xorl	%r8d, %r8d
	movq	%rbx, %rdi
	movsd	%xmm1, (%rsp)
	call	assign_forces_ad
	call	MPI_Wtime
	movsd	(%rsp), %xmm1
	movl	$.LC8, %edi
	movl	$1, %eax
	subsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movsd	%xmm1, (%rsp)
	call	printf
	movl	$.LC9, %edi
	call	puts
	movsd	16(%rsp), %xmm2
	movl	$.LC11, %edi
	movsd	.LC10(%rip), %xmm3
	movl	$1, %eax
	movapd	%xmm2, %xmm0
	movapd	%xmm3, 32(%rsp)
	xorpd	%xmm3, %xmm0
	call	printf
	movsd	(%rsp), %xmm1
	movl	$.LC12, %edi
	movsd	56(%rsp), %xmm0
	movl	$1, %eax
	movapd	32(%rsp), %xmm3
	addsd	%xmm1, %xmm0
	xorpd	%xmm3, %xmm0
	call	printf
	movq	%rbp, %rdi
	call	Free_data
	addq	$112, %rsp
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE19:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	1076101120
	.align 8
.LC1:
	.long	0
	.long	1072693248
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC10:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.ident	"GCC: (SUSE Linux) 4.6.2"
	.section	.comment.SUSE.OPTs,"MS",@progbits,1
	.string	"Ospwg"
	.section	.note.GNU-stack,"",@progbits
