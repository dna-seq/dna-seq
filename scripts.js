/* let els = document.querySelectorAll('.sections')
  els = Array.from(els)


const observer = new IntersectionObserver(entries => {
    entries.forEach(entry => {
      const square = entry.target.querySelector('.test');
  
      if (entry.isIntersecting) {
        square.classList.add('sections-animation');
        return; // if we added the class, exit the function
      }
  
      // We're not intersecting, so remove the class!
      square.classList.remove('sections-animation');
    });
  });
  
  observer.observe(document.querySelector('.wrapper')); */

/*  function scrollTrigger(selector){
    let els = document.querySelectorAll(selector)
    els = Array.from(els)
    els.forEach(el => {
      addObserver(el)
    })
  }

  function addObserver(el, options){
    let observer = new IntersectionObserver((entries, observer) => {
      entries.forEach(entry => {
        if(entry.isIntersecting){
          if(options.cb) {
            // If we've passed a callback function, we'll call it
            options.cb(el)
          } else{
            // If we haven't, we'll just add the active class
            entry.target.classList.add('sections-animation')
          }
          observer.unobserve(entry.target)
        }
      })
    }, options)
    observer.observe(el)
  }

  scrollTrigger('.sections1', {
    rootMargin: '-200px'
  })*/

const observer = new IntersectionObserver((entries) => {
  entries.forEach((entry) => {
    console.log(entry)
    if (entry.isIntersecting) {
      entry.target.classList.add('.sections-animation');
    } else {
      entry.target.classList.remove('.sections-animation');
    }
  });
});

const hiddenElements = document.querySelectorAll('.sections');
hiddenElements.forEach((el) => observer.observe(el));
