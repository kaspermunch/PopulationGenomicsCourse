---
title: General page
summary: A brief description of my document.
date: 2000-01-01
hide:
  - navigation
---

<!--
# Put above to hide navigation (left), toc (right) or footer (bottom)

hide:
  - navigation 
  - toc
  - footer 

# You should hide the navigation if there are no subsections
# You should hide the Table of Contents if there are no important titles
-->

# General Page

!!! note "Section Overview"

    &#128368; **Time Estimation:** X minutes  

    &#128172; **Learning Objectives:**    
        1. First item  
        2. Second item  

    
Write your introduction to the page here.

<hr>

## Useful functions

You should continue to write your markdown document as normal. But here are some useful functions. 
You can find more in the Reference guide listed on the tab above. **Be sure to delete the guide 
and the tips below when finished with them.**

### Code

Text will be highlighted appropriately when you include language abbreviation:

```py
import tensorflow as tf
```

### Admonitions

The admonitions used for course/section overview and requirements should be consistent, though you 
can use any other admonitions freely.

Examples (both drop-down and not):

<!-- !!! = no drop-down -->

!!! quote
    Here is a quote

<!-- ??? = drop-down -->

??? question "What is the smallest country in the world?"
    A: Vatican City


### Footnotes

We can include footnotes like this one[^1].


### LaTeX

You write an equation as normal:

$$
\operatorname{ker} f=\{g\in G:f(g)=e_{H}\}{\mbox{.}}
$$


### Images 

Format is similar to links, but include an exclamation mark before: 

![Image title](https://dummyimage.com/600x400/eee/aaa)

You can link to a URL or to somewhere locally.

<!-- Footnote content -->

[^1]: Remember to eat your vegetables. 