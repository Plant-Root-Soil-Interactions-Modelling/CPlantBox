/*!
* Start Bootstrap - Agency v7.0.12 (https://startbootstrap.com/theme/agency)
* Copyright 2013-2023 Start Bootstrap
* Licensed under MIT (https://github.com/StartBootstrap/startbootstrap-agency/blob/master/LICENSE)
*/
//
// Scripts
// 

window.addEventListener('DOMContentLoaded', event => {

    // Navbar shrink function
    var navbarShrink = function () {
        const navbarCollapsible = document.body.querySelector('#mainNav');
        if (!navbarCollapsible) {
            return;
        }
        if (window.scrollY === 0) {
            navbarCollapsible.classList.remove('navbar-shrink')
        } else {
            navbarCollapsible.classList.add('navbar-shrink')
        }

    };

    // Shrink the navbar 
    navbarShrink();

    // Shrink the navbar when page is scrolled
    document.addEventListener('scroll', navbarShrink);

    //  Activate Bootstrap scrollspy on the main nav element
    const mainNav = document.body.querySelector('#mainNav');
    if (mainNav) {
        new bootstrap.ScrollSpy(document.body, {
            target: '#mainNav',
            rootMargin: '0px 0px -40%',
        });
    };

    // Collapse responsive navbar when toggler is visible
    const navbarToggler = document.body.querySelector('.navbar-toggler');
    const responsiveNavItems = [].slice.call(
        document.querySelectorAll('#navbarResponsive .nav-link')
    );
    responsiveNavItems.map(function (responsiveNavItem) {
        responsiveNavItem.addEventListener('click', () => {
            if (window.getComputedStyle(navbarToggler).display !== 'none') {
                navbarToggler.click();
            }
        });
    });

    // Gallery filters and image preview modal
    const galleryFilterButtons = [].slice.call(
        document.querySelectorAll('[data-gallery-filter]')
    );
    const galleryItems = [].slice.call(
        document.querySelectorAll('.cpb-gallery-item')
    );

    galleryFilterButtons.forEach((button) => {
        button.addEventListener('click', () => {
            const filter = button.dataset.galleryFilter;

            galleryFilterButtons.forEach((filterButton) => {
                const isActive = filterButton === button;
                filterButton.classList.toggle('active', isActive);
                filterButton.classList.toggle('btn-primary', isActive);
                filterButton.classList.toggle('btn-outline-primary', !isActive);
                filterButton.setAttribute('aria-pressed', isActive.toString());
            });

            galleryItems.forEach((item) => {
                const matchesCategory = filter === 'all' || item.dataset.galleryCategory === filter;
                const matchesHighlight = filter === 'highlight' && item.dataset.galleryFeatured === 'true';
                item.classList.toggle('is-hidden', !(matchesCategory || matchesHighlight));
            });
        });
    });

    const galleryModalElement = document.getElementById('galleryModal');
    if (galleryModalElement) {
        const galleryModal = new bootstrap.Modal(galleryModalElement);
        const galleryModalImage = galleryModalElement.querySelector('.cpb-gallery-modal-image');
        const galleryModalTitle = galleryModalElement.querySelector('#galleryModalLabel');
        const galleryModalPrev = galleryModalElement.querySelector('.cpb-gallery-modal-prev');
        const galleryModalNext = galleryModalElement.querySelector('.cpb-gallery-modal-next');
        let currentGalleryIndex = 0;

        const getVisibleGalleryLinks = () => [].slice.call(
            document.querySelectorAll('.cpb-gallery-item:not(.is-hidden) .cpb-gallery-link')
        );

        const showGalleryImage = (index) => {
            const visibleGalleryLinks = getVisibleGalleryLinks();
            if (visibleGalleryLinks.length === 0) {
                return;
            }

            currentGalleryIndex = (index + visibleGalleryLinks.length) % visibleGalleryLinks.length;
            const galleryLink = visibleGalleryLinks[currentGalleryIndex];

            galleryModalImage.src = galleryLink.href;
            galleryModalImage.alt = galleryLink.querySelector('img').alt;
            galleryModalTitle.textContent = galleryLink.dataset.galleryTitle;
        };

        const stepGalleryImage = (step) => {
            showGalleryImage(currentGalleryIndex + step);
        };

        document.querySelectorAll('.cpb-gallery-link').forEach((galleryLink) => {
            galleryLink.addEventListener('click', (event) => {
                event.preventDefault();

                currentGalleryIndex = getVisibleGalleryLinks().indexOf(galleryLink);
                showGalleryImage(currentGalleryIndex);
                galleryModal.show();
            });
        });

        galleryModalPrev.addEventListener('click', () => {
            stepGalleryImage(-1);
        });

        galleryModalNext.addEventListener('click', () => {
            stepGalleryImage(1);
        });

        galleryModalElement.addEventListener('keydown', (event) => {
            if (event.key === 'ArrowLeft') {
                event.preventDefault();
                stepGalleryImage(-1);
            }

            if (event.key === 'ArrowRight') {
                event.preventDefault();
                stepGalleryImage(1);
            }
        });
    }

});
